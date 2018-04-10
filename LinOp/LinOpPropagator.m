classdef LinOpPropagator <  LinOp
    %% LinOpPropagator: propagation operator
    %  Propagation linear operator that model the scalar propagation of coherent light through an homogene medium
    %
    % Example
    % Obj = LinOpPropagator(sz,lambda, n0, z,dxy, theta, illu,type,NA)
    % Fresnel transform operator with
    %  :param sz:       size of the input
    %  :param lambda:   wavelenght [m]
    %  :param n0:       refractive index of the medium
    %  :param z:        depth of propagation  [m]
    %  :param dxy:      pixel size            [m]
    %  :param theta:    incidence angle along x and y. It is a vector [2  Nt] where Nt is the number of incidence angles.
    %  :param illu:     the illumination wave in the sample plane. It can be as scalar for uniform illumination, a vector [1 Nt] for a different but uniform illumination at each incidence angle or a (complex) vector [sz Nt] for a random illumination
    %  :param type:     model of the diffraction pattern. It can be either:
    %                        * 'Fresnel' for Fresnel method (in Fourier domain)
    %                        * 'AS'      for Angular Spectrum method
    %                        * 'FeitFleck' for the Feit and Fleck model   (M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.)
    %                        * 'Pupil'   for propagation through an objective with from the focal plane to the detector plane. In that case the next argument is :param NA: the numerical apperture.
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %     Copyright (C) 2018 F. Soulez ferreol.soulez@epfl.ch
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
        FRESNEL = 0
    end
    properties (SetAccess = protected,GetAccess = public)
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        Nt      % number of angle
        ephi       % Fresnel function
        sizeinC;
        sizeoutC;
        unifInten = true;
        modulus;
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
        function this = LinOpPropagator(sz,lambda, n0, z,dxy,theta, illu,varargin)
            
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
            this.sizeinC = this.sizein;
            this.sizeoutC = this.sizeout;
            
            this.Nx = sz(1);
            this.Ny = sz(2);
            
            
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
                end
            end
            
            
            
            this.update();
            
            addlistener(this,{'theta','z','lambda','type','n0','NA','dxy','illu'},'PostSet',@this.update);
            
        end
    end
    methods  (Access = protected)
        function y = apply_(this,x)
            y = zeros(this.sizeoutC);
            
            for nt = 1:this.Nt
                fx = fft2(this.illu(:,:,nt).*x);
                y(:,:,nt) = ifft2( this.ephi(:,:,nt) .*  fx);
            end
        end
        function y = applyAdjoint_(this,x)
            y = zeros(this.sizeinC);
            for nt = 1:this.Nt
                y = y+ conj(this.illu(:,:,nt)).*  ifft2( conj(this.ephi(:,:,nt)) .*  fft2(x(:,:,nt)));
            end
            
        end
        function y = applyHtH_(this, x)
            
                 y = x.* sum(abs(this.illu).^2,3);
        end        
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            M=LinOpDiag(this.sizein,sum(abs(this.illu).^2,3));
        end
    end
    methods (Access = private)
        function update(this,~,~)
            
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            dxy_ = this.dxy;
            theta_ = this.theta.*ones(2,this.Nt);
            n0_ = this.n0;
            lambda_ = this.lambda.*ones(1,this.Nt);
            z_ = this.z.*ones(1,this.Nt);
            
            
            this.unifInten =  true;
            
            % propagation kernel
            ephi_ = zeros(this.sizeoutC);
            
            if this.type==this.PUPIL
                mod = ones([this.sizein, this.Nt]);
            else
                mod = zeros([1,1,this.Nt]);
            end
            for nt = 1:this.Nt
                %  frequency grid
                v = 1./(Nx_ * dxy_) *( [0:ceil( Nx_/2)-1, -floor( Nx_/2):-1]' -   dxy_ *   sin(theta_(1,nt))/lambda_(nt).*Nx_);
                u = 1./(Ny_ * dxy_) * ([0:ceil( Ny_/2)-1, -floor( Ny_/2):-1] -    dxy_ *   sin(theta_(2,nt))/lambda_(nt).*Ny_);
                
                
                [mu,mv] =  meshgrid(  u.^2,  v.^2);
                Mesh = mv + mu;
                
                if this.type~=this.PUPIL
                    if (max(Mesh(:))>(n0_/ lambda_(nt))^2)
                        if isequal(size(mod),[1,1,this.Nt])
                            mod = ones([this.sizein this.Nt]);
                        end
                        mod(:,:,nt) =  (Mesh<= (n0_/ lambda_(nt))^2);
                        Mesh(~mod(:,:,nt))=0.;
                        this.unifInten =  false;
                    else
                        mod(:,:,nt) = 1;
                    end
                end
                
                switch  this.type
                    case  this.PUPIL % pupil
                        mod(:,:,nt) = (Mesh<=(this.NA/ lambda_(nt))^2);
                        ephi_(:,:,nt) = mod(:,:,nt);
                        this.unifInten =  false;
                    case  this.AS % Angular spectrum
                        ephi_(:,:,nt) =  mod(:,:,nt).*exp(-2i* pi *  z_(nt).* sqrt((  n0_/ lambda_(nt))^2- Mesh));
                    case  this.FEITFLECK
                        ephi_(:,:,nt)=  mod(:,:,nt).*exp(-2i* pi * z_(nt).*lambda_(nt) / n0_ * Mesh ./ real(1 + sqrt(1 - (lambda_(nt)/n0_)^2 *Mesh)));
                    otherwise
                        %  Fresnel function
                        ephi_(:,:,nt) =mod(:,:,nt).* exp(-1i* pi *  z_(nt).* lambda_(nt) / n0_ .*Mesh);
                end
            end
            this.ephi = ephi_;
            this.modulus = mod.^2;
            
            % Illumination
            
            if isscalar(this.illu) || (numel(this.illu)==this.Nt)
                this.illu =  ones([1,1,this.Nt]).*this.illu;
            elseif this.unifInten
                for nt = 1:this.Nt
                    tmp = abs(this.illu(:,:,nt));
                    if( length(unique( tmp(:)))~=1)
                        this.unifInten =  false;
                        break;
                    end
                end
            end
            
            
        end
    end
end
