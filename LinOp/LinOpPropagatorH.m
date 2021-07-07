classdef LinOpPropagatorH <  LinOp
    %% LinOpPropagatorH : propagation operator in Fourier
    %  Linear operator that models the scalar propagation of coherent light through an homogene medium
    %  It is identical to the :class:`LinOpPropagator``excepted that its input is already in Fourier domain bypassing some extra FFT calls
    %
    %  :param sz:       size of the input
    %  :param lambda:   wavelenght [m] can be a scalar or
    %                   a vector of size [2  Nt] where Nt is the number image.
    %  :param n0:       refractive index of the medium can be a scalar or
    %                   a vector of size [2  Nt] where Nt is the number image.
    %  :param z:        depth of propagation vector  [m] can be a scalar or
    %                   a vector of size [2  Nt] where Nt is the number image.
    %  :param dxy:      pixel size            [m]
    %  :param theta:    incidence angle along x and y. It is a vector [2  Nt] where Nt is the number of incidence angles.
    %  :param illu:     the illumination wave in the sample plane. It can be either a scalar for uniform illumination or a vector [1 Nt] for a different but uniform illumination at each incidence angle
    %  :param fourierWeight:   array of [sizeout] for an extra filtering in Fourier to model the finite coherence,
    %  :param type:     model of the diffraction pattern. It can be either:
    %                        * 'Fresnel' for Fresnel method (in Fourier domain)
    %                        * 'AS'      for Angular Spectrum method
    %                        * 'BLAS'    for Band Limited Angular Spectrum (Matsushima, K., & Shimobaba, T. (2009).   Band-limited angular spectrum method for numerical simulation of free-space propagation in far and near fields. Optics express, 17(22), 19662-19673.)
    %                        * 'BLAS2'   idem Yu, X., Xiahui, T., Yingxiong, Q., Hao, P., & Wei, W. (2012). Band-limited angular spectrum numerical propagation method with selective scaling of observation window size and sample number. JOSA A, 29(11), 2415-2420.
    %                        * 'FeitFleck' for the Feit and Fleck model   (M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.)
    %                        * 'Pupil'   for propagation through an objective with from the focal plane to the detector plane. In that case the next argument is :param NA: the numerical apperture.
    %  :param PhaseCorrect if the keyword 'PhaseCorrect' is set, the output will be modulated by exp(1I sin(theta)/lambda) due to the tilted incidence. This therm is usually not needed as it vanished from the intensity
    %  :param oversample if the keyword 'oversample' is set, the propagation is oversampled by a factor 2 to prevent aliasing in subsequent intensity estimation. The output size will be twice the input size.
    %
    %
    %  All attributes of parent class :class:`LinOp` are inherited.
    %
    %  **Example**  H = LinOpPropagatorH(sz,lambda, n0, z,dxy, theta, illu,type,NA,oversample)
    %
    % See also :class:`LinOpPropagator`, :class:`LinOp`, :class:`Map`
    %
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
    
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (Constant=true)
        PUPIL = 3   % Fresnel propagation through an objective of pupil of radius R
        FEITFLECK = 2
        AS = 1 % Use Angular Spectrum method
        BLAS = 4; % Band limited AS (
        BLAS2 = 5;
        FRESNEL = 0; % Fresnel propagation model
    end
    properties (SetAccess = protected,GetAccess = public)
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        Nt      % number of angle
        ephi       % Fresnel function
        unifFourier =  true % true if the kernel has a unit modulus
        fourierWeight =1; %  weight in Fourier space (eg attenuation given by the finite coherence)
        N       % number of pixel
        scale;  % global illumination scale per hologram
        PhRampV;
        PhRampU;
        oversample= false;
        sizeoutno % sizeout without oversampling
    end
    
    properties (SetObservable)
        type = 0;
        lambda  % wavelenght            [m]
        n0      % refractive index of the medium
        z       % depth of propagation  [m]
        dxy     % pixel size            [m]
        theta   % tilt angle
        NA =  1;% numerical aperture (if needed)
        PhaseCorrect= false
    end
    methods
        function this = LinOpPropagatorH(sz,lambda, n0, z,dxy, theta,illu,fourierWeight,varargin)
            
            this.name ='LinOpPropagatorH';
            
            this.lambda = lambda;
            
            assert(isPositiveScalar(n0),'The refractive index n0 should be a positive scalar');
            this.n0 = n0;
            
            this.z = z;
            
            assert(isPositiveScalar(dxy),'The pixel size dxy should be a positive scalar');
            this.dxy = dxy;
            
            assert(issize(sz) && (length(sz)==2),'The input size sz should be a conformable  to size(2D) ');
            
            this.theta = theta;
            this.Nt = max([numel(theta)/2,numel(lambda),numel(z)]);
            this.sizein = sz ;
            this.sizeout = [sz this.Nt];
            this.sizeoutno = [sz this.Nt];
            
            this.Nx = sz(1);
            this.Ny = sz(2);
            this.N = prod(sz);
            this.scale=ones_(1,1 , this.Nt);
            
            this.fourierWeight = fourierWeight;
            for c=1:length(varargin)
                switch varargin{c}
                    case('Pupil')
                        this.type = this.PUPIL;
                        this.NA = varargin{c+1};
                    case('FeitFleck')
                        this.type = this.FEITFLECK;
                    case('AS')
                        this.type = this.AS;
                    case('BLAS')
                        this.type = this.BLAS;
                    case('BLAS2')
                        this.type = this.BLAS2;
                    case('Fresnel')
                        this.type = this.FRESNEL;
                    case('PhaseCorrect')
                        this.PhaseCorrect=true;
                    case('oversample')
                        this.oversample=true;
                        this.sizeout = [2.*sz this.Nt];
                    otherwise
                        error('Type of propagation must be defined');
                end
            end
            
            % Illumination
            if isscalar(illu) || (numel(illu)==this.Nt)
                this.scale = this.scale.* illu(:);
            else
                error('illumination must be uniform');
            end
            
            this.update();
            
            addlistener(this,{'theta','z','lambda','type','n0','NA','dxy','PhaseCorrect'},'PostSet',@this.update);
        end
    end
    methods  (Access = protected)
        function y = apply_(this,x)
            y = zeros_(this.sizeout);           
            for nt = 1:this.Nt
                fx =this.ephi(:,:,nt) .*x;
                if this.oversample
                    fx =2.*ifft2(ifftshift( padarray( fftshift(fx),[this.Nx/2, this.Ny/2] )));
                else
                    fx =ifft2(fx);
                end
                if this.PhaseCorrect
                    y(:,:,nt) = this.scale(nt).*this.PhRampV(:,nt)'.*this.PhRampU(:,nt).*  fx;
                else
                    y(:,:,nt) = this.scale(nt).*fx;
                end
            end
            
        end
        function y = applyAdjoint_(this,x)
            y = zeros_(this.sizein);
            
            for nt = 1:this.Nt
                fx = fft2(x(:,:,nt));
                if this.oversample
                    fx = circshift(fx,[this.Nx/2, this.Ny/2] );
                    fx = 0.5*ifftshift( fx(1:this.Nx, 1:this.Ny));
                end
                y = y+ this.scale(nt)./(this.N) .*    conj(this.ephi(:,:,nt)) .* fx;
            end
            
        end
        function y = applyHtH_(this,x)
            m = real(this.ephi).^2 + imag(this.ephi).^2;
            sc2 = abs(this.scale).^2;
            y = sum(sc2.*m,3)./(this.N).*x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            m = real(this.ephi).^2 + imag(this.ephi).^2;
            sc2 = abs(this.scale).^2;
            M=LinOpDiag(this.sizein,sum(sc2.*m,3)./(this.N) );
            
            
        end
    end
    
    methods (Access = private)
        function update(this,~,~)
            
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            dxy_ = this.dxy;
            
            theta_ = this.theta.*ones_(2,this.Nt);
            n0_ = this.n0;
            lambda_ = this.lambda.*ones_(1,this.Nt);
            z_ = this.z.*ones_(1,this.Nt);
            
            
            % Fourier Weight
            if isscalar(this.fourierWeight) || numel(this.fourierWeight)==this.Nt
                fourierWeight_ = ones_(1,1, this.Nt);
                this.scale = this.scale.*this.fourierWeight(:);
                this.fourierWeight = 1;
            elseif cmpSize(this.sizeoutno, size(this.fourierWeight))
                fourierWeight_ = this.fourierWeight;
                for nt = 1:this.Nt
                    tmp = abs(this.fourierWeight(:,:,nt));
                    if( length(unique( tmp(:)))~=1)
                        this.unifFourier =  false;
                        break;
                    end
                end
            else
                error('fourierWeight does not have the right size');
            end
            
            
            
            
            
            % propagation kernel
            ephi_ = zeros_(this.sizeoutno);
            if this.PhaseCorrect
                if this.oversample
                    this.PhRampV = zeros_(2.*Nx_,this.Nt);
                    this.PhRampU = zeros_(2.*Ny_,this.Nt);
                else
                    this.PhRampV = zeros_(Nx_,this.Nt);
                    this.PhRampU = zeros_(Ny_,this.Nt);
                end
            end
            
            
            L = min(Nx_,Ny_)* dxy_/2;
            
            for nt = 1:this.Nt
                %  frequency grid
                v = 1./(Nx_ *  dxy_) *( [0:ceil( Nx_/2)-1, -floor( Nx_/2):-1]' -    dxy_ *   sin( theta_(1,nt)) *n0_/ lambda_(nt).*Nx_);
                u = 1./(Ny_ *  dxy_) * ([0:ceil( Ny_/2)-1, -floor( Ny_/2):-1] -     dxy_ *   sin( theta_(2,nt)) *n0_/ lambda_(nt).*Ny_);
                
                
                [mu,mv] =  meshgrid(  u.^2,  v.^2);
                Mesh = mv + mu;
                
                if this.type==this.BLAS
                    Q = (L/sqrt(L^2+z_(nt)^2) *n0_/ lambda_(nt))^2;
                    if (max(Mesh(:))>Q)
                        mod  =  fourierWeight_(:,:,nt).*(mv<= Q*(1-mu* ((lambda_(nt)/ n0_)^2))).*(mu<= Q*(1-mv* ((lambda_(nt)/ n0_)^2)));
                        Mesh(~mod )=0.;
                        this.unifFourier =  false;
                    else
                        mod = fourierWeight_(:,:,nt);
                    end
                else
                    if this.type==this.BLAS2
                        Q = (L/sqrt(L^2+z_(nt)^2) *n0_/ lambda_(nt))^2;
                    else
                        Q=(n0_/ lambda_(nt))^2;
                    end
                    if (max(Mesh(:))>Q)
                        mod  =  fourierWeight_(:,:,nt).*(Mesh<= Q);
                        Mesh(~mod )=0.;
                        this.unifFourier =  false;
                    else
                        mod = fourierWeight_(:,:,nt);
                    end
                end
                clear mv mu;
                switch  this.type
                    case  this.PUPIL % pupil
                        ephi_(:,:,nt) = mod.*(Mesh<=(this.NA/ lambda_(nt))^2);
                        this.unifFourier =  false;
                    case  {this.AS,this.BLAS,this.BLAS2} % Angular spectrum
                        ephi_(:,:,nt) =  mod.*exp(2i* pi *  z_(nt).* sqrt((  n0_/ lambda_(nt))^2- Mesh));
                    case  this.FEITFLECK
                        ephi_(:,:,nt)=  mod.*exp(-2i* pi * z_(nt).*lambda_(nt) / n0_ * Mesh ./ real(1 + sqrt(1 - (lambda_(nt)/n0_)^2 *Mesh)));
                    otherwise
                        %  Fresnel function
                        ephi_(:,:,nt) =mod.* exp(-1i* pi *  z_(nt).* lambda_(nt) / n0_ .*Mesh);
                end
                
                
                
                this.norm=max(max(this.scale(nt).*abs(ephi_(:))),this.norm);
                
                
                if this.PhaseCorrect
                    if this.oversample
                        this.PhRampV(:,nt) = exp( (2.i * pi *dxy_ *   sin(theta_(1,nt)) *n0_/lambda_(nt)) .* (1:2*Nx_));
                        this.PhRampU(:,nt) = exp( (2.i * pi *dxy_ *   sin(theta_(2,nt)) *n0_/lambda_(nt)) .* (1:2*Ny_));
                    else
                        this.PhRampV(:,nt) = exp( (2.i * pi *dxy_ *   sin(theta_(1,nt)) *n0_/lambda_(nt)) .* (1:Nx_));
                        this.PhRampU(:,nt) = exp( (2.i * pi *dxy_ *   sin(theta_(2,nt)) *n0_/lambda_(nt)) .* (1:Ny_));
                    end
                end
            end
            this.ephi = ephi_;
        end
    end
end
