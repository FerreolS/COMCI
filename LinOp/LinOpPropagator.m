classdef LinOpPropagator <  LinOp
    %% LinOpPropagator: propagation operator
    %  Propagation linear operator that model the scalar propagation of coherent light through an homogene medium
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
    %  :param illu:     the illumination wave in the sample plane. It can be as scalar for uniform illumination, a vector [1 Nt] for a different but uniform illumination at each incidence angle or a (complex) vector [sz Nt] for a random illumination
    %  :param fourierWeight:   array of [sizeout] for an extra filtering in Fourier to model the finite coherence,
    %  :param type:     model of the diffraction pattern. It can be either:
    %                        * 'Fresnel' for Fresnel method (in Fourier domain)
    %                        * 'AS'      for Angular Spectrum method
    %                        * 'BLAS'    for Band Limited Angular Spectrum (Matsushima, K., & Shimobaba, T. (2009).   Band-limited angular spectrum method for numerical simulation of free-space propagation in far and near fields. Optics express, 17(22), 19662-19673.)
    %                        * 'BLAS2'   idem Yu, X., Xiahui, T., Yingxiong, Q., Hao, P., & Wei, W. (2012). Band-limited angular spectrum numerical propagation method with selective scaling of observation window size and sample number. JOSA A, 29(11), 2415-2420.
    %                        * 'FeitFleck' for the Feit and Fleck model   (M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.)
    %                        * 'Pupil'   for propagation through an objective with from the focal plane to the detector plane. In that case the next argument is :param NA: the numerical apperture.
    %
    %  All attributes of parent class :class:`LinOp` are inherited.
    %
    %  **Example**  H = LinOpPropagator(sz,lambda, n0, z,dxy, theta, illu,type,NA)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %     Copyright (C) 2018 F. Soulez ferreol.soulez@univ-lyon1.fr
    %
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
        FRESNEL = 0
        
    end
    properties (SetAccess = protected,GetAccess = public)
        Nx      % number of pixels along X            
        Ny      % number of pixels along Y
        Nt      % number of angle
        ephi       % Fresnel function
        unifIllu = true; % true if the illumination intensity is uniform
        unifFourier =  true % true if the kernel has a unit modulus
        fourierWeight =1; %  weight in Fourier space (eg attenuation given by the finite coherence)
        N       % number of pixel
        scale;  % global illumination scale per hologram
        
    end
    properties (SetObservable)
        type = 0;
        illu;
        lambda  % wavelenght            [m]
        n0      % refractive index of the medium
        z       % depth of propagation  [m]
        dxy     % pixel size            [m]
        theta   % tilt angle
        NA =  1;% numerical aperture (if needed)
    end
    methods
        function this = LinOpPropagator(sz,lambda, n0, z,dxy,theta, illu,fourierWeight,varargin)
            
            this.name ='LinOpPropagator';
            
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
            
            this.Nx = sz(1);
            this.Ny = sz(2);
            this.N = prod(sz);
            this.scale=ones(1,1 , this.Nt);
            
            this.fourierWeight = fourierWeight;
            this.illu = illu;
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
                    otherwise
                        error('Type of propagation must be defined');
                end
            end
            
            
            
            this.update();
            
            
            addlistener(this,{'theta','z','lambda','type','n0','NA','dxy','illu'},'PostSet',@this.update);
            
        end
    end
    methods  (Access = protected)
        function y = apply_(this,x)
            y = zeros(this.sizeout);
            
            for nt = 1:this.Nt
                fx = fft2(this.illu(:,:,nt).*x);
                y(:,:,nt) = this.scale(nt).*ifft2( this.ephi(:,:,nt) .*  fx);
            end
        end
        function y = applyAdjoint_(this,x)
            y = zeros(this.sizein);
            for nt = 1:this.Nt
                y = y+ this.scale(nt).*conj(this.illu(:,:,nt)).*  ifft2( conj(this.ephi(:,:,nt)) .*  fft2(x(:,:,nt)));
            end
            
        end
        function y = applyHtH_(this, x)
            if this.unifFourier
                y = x.* sum(abs(this.illu).^2,3).*mean(this.scale(:).^2);
            elseif this.unifIllu
                m = real(this.ephi).^2 + imag(this.ephi).^2;
                y =ifftn(sum(m,3).*mean(this.scale(:).^2).*fftn(x));
            else
                y= applyHtH_@LinOp(this,x);
            end
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unifFourier
                M=LinOpDiag(this.sizein,sum(abs(this.illu).^2,3) .*mean(this.scale(:).^2));
            elseif this.unifIllu
                m = real(this.ephi).^2 + imag(this.ephi).^2;
                M=LinOpConv(sum(m,3).*mean(this.scale(:).^2),false );
            else
                M=     makeHtH_@LinOp(this);
            end
        end
    end
    methods (Access = private)
        function update(this,~,~)
            
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            dxy_ = this.dxy;
            if this.Nt==1 && size(this.theta,1)==1
                this.theta = this.theta';
            end
            theta_ = this.theta.*ones(2,this.Nt);
            n0_ = this.n0;
            lambda_ = this.lambda.*ones(1,this.Nt);
            z_ = this.z.*ones(1,this.Nt);
            
            
            
            % Fourier Weight
            if isscalar(this.fourierWeight) || numel(this.fourierWeight)==this.Nt
                fourierWeight_ = ones(1,1, this.Nt);
                this.scale = this.scale.*this.fourierWeight(:);
                this.fourierWeight = 1;
            elseif cmpSize(this.sizeout, size(this.fourierWeight))
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
            
            
            % Illumination
            this.unifIllu =  true;
            if isscalar(this.illu) || (numel(this.illu)==this.Nt)
                this.scale = this.scale.* this.illu(:);
                this.illu =  ones([1,1,this.Nt]);
            elseif cmpSize(this.sizeout, size(this.illu))
                for nt = 1:this.Nt
                    tmp = abs(this.illu(:,:,nt));
                    if( length(unique( tmp(:)))~=1)
                        this.unifIllu =  false;
                        break;
                    else
                        this.scale(nt) = tmp;
                        this.illu(:,:,nt) = this.illu(:,:,nt) ./tmp;
                    end
                end
            else
                error('illu does not have the right size');
            end
            
            
            
            % propagation kernel
            ephi_ = zeros(this.sizeout);
            
            L = Nx_* dxy_/2;
            for nt = 1:this.Nt
                %  frequency grid
                v = 1./(Nx_ * dxy_) *( [0:ceil( Nx_/2)-1, -floor( Nx_/2):-1]' -   dxy_ *   sin(theta_(1,nt))/lambda_(nt).*Nx_);
                u = 1./(Ny_ * dxy_) * ([0:ceil( Ny_/2)-1, -floor( Ny_/2):-1] -    dxy_ *   sin(theta_(2,nt))/lambda_(nt).*Ny_);
                
                
                [mu,mv] =  meshgrid(  u.^2,  v.^2);
                Mesh = mv + mu;
                if this.type==this.BLAS
                    Q = (L/sqrt(L^2+z_(nt)^2) *n0_/ lambda_(nt))^2;
                    if (max(Mesh(:))>Q)
                        mod  =  fourierWeight_(:,:,nt).*(mv<= Q*(1-mu* (lambda_(nt)^2))).*(mu<= Q*(1-mv* (lambda_(nt)^2)));
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
                    case  {this.AS,this.BLAS} % Angular spectrum
                        ephi_(:,:,nt) =  mod.*exp(-2i* pi *  z_(nt).* sqrt((  n0_/ lambda_(nt))^2- Mesh));
                    case  this.FEITFLECK
                        ephi_(:,:,nt)=  mod.*exp(-2i* pi * z_(nt).*lambda_(nt) / n0_ * Mesh ./ real(1 + sqrt(1 - (lambda_(nt)/n0_)^2 *Mesh)));
                    otherwise
                        %  Fresnel function
                        ephi_(:,:,nt) =mod.* exp(-1i* pi *  z_(nt).* lambda_(nt) / n0_ .*Mesh);
                end
            end
            this.ephi = ephi_;
            
            
            this.norm=max(this.scale(nt).*abs(this.ephi(:)));
        end
    end
end
