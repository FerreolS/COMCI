classdef LinOpPropagatorH <  LinOp
    %% LinOpPropagatorH : Stack of tilted propagator operator in Fourier
    %
    %  Matlab Linear Operator Library
    %    % Fresnel transform operator for a stack of T holograms with
    %    :param lambda:   array of T wavelenghts in   [m]
    %    :param     n0:   refractive index of the medium
    %    :param      z:   array of T depth of propagation  [m]
    %    :param    dxy:   pixel size (width=lenght =dxy)   [m]
    %    :param  theta:   array of [2 T] illumination angle in radian
    %    :param   illu:   array of T mean amplitude of the illumination per hologram 
    %    :param   type:   type of the proagation model can be:
    %                      'Fresnel' (default)
    %                      'FeitFleck' use the Feit and Fleck model (M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.)
    %                      'AS' for angular spectrum 
    %                      'Pupil' for propagation through an objective
    %         if the option 'Pupil' is set then the extra :param NA: for numerical apperture is needed
    %
    % Example
    % Obj =  LinOpStackOfTiltedPropagatorH(sz,lambda, n0, z,dxy, theta, illu,'AS')
    %
    % Obj =  LinOpStackOfTiltedPropagatorH(sz,lambda, n0, z,dxy, theta, illu,'Pupil',NA)
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
        FRESNEL = 0 % Fresnel propagation model
        PUPIL = 3   % Fresnel propagation through an objective of pupil of numerical apperture NA
        FEITFLECK = 2 % use the Feit and Fleck model of propagation instead:
        % M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
        AS = 1 % Use Angular Spectrum method
    end
    properties (SetAccess = protected,GetAccess = public)
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        Nt      % number of angle
        ephi       % Fresnel function
        sizeinC;
        sizeoutC;
        
    end
    
    properties (SetObservable)
        type = 0;
        lambda  % wavelenght            [m]
        n0      % refractive index of the medium
        z       % depth of propagation  [m]
        dxy     % pixel size            [m]
        theta   % tilt angle
        NA =  1;% numerical aperture (if needed)
        illu=1; %illumination
        N       % number of pixel
    end
    methods
        function this = LinOpPropagatorH(sz,lambda, n0, z,dxy, theta, illu,varargin)
            
            this.name ='LinOpPropagatorH';
            
            this.lambda = lambda;
            
            assert(isPositiveScalar(n0),'The refractive index n0 should be a positive scalar');
            this.n0 = n0;
            
            this.z = z;
            
            assert(isPositiveScalar(dxy),'The pixel size dxy should be a positive scalar');
            this.dxy = dxy;
            
            assert(issize(sz) && (length(sz)==2),'The input size sz should be a conformable  to size(2D) ');
            
            this.theta = theta;
            this.Nt = numel(theta)/2;
            this.sizein = sz ;
            this.sizeout = [sz this.Nt];
            this.sizeinC = this.sizein;
            this.sizeoutC = this.sizeout;
            
            this.Nx = sz(1);
            this.Ny = sz(2);
            
            
            this.N = prod(sz);
            
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
            
            addlistener(this,{'theta','z','lambda','type','n0','NA','dxy'},'PostSet',@this.update);
        end
    end
    methods  (Access = protected)
        function y = apply_(this,x)
            y = zeros(this.sizeoutC);
            for nt = 1:this.Nt
                y(:,:,nt) = ifft2( this.ephi(:,:,nt) .*  x);
            end
        end
        function y = applyAdjoint_(this,x)
            y = zeros(this.sizeinC);
            for nt = 1:this.Nt
                y = y+    conj(this.ephi(:,:,nt)) .*  fft2(x(:,:,nt));
            end
            y = y.* 1./(this.N) ;
        end
        function y = applyHtH_(this,x)
            re = real(this.ephi).^2 + imag(this.ephi).^2;
            y = sum(re,3)./(this.N).*x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            re = real(this.ephi).^2 + imag(this.ephi).^2;
            
            M=LinOpDiag(this.sizein,sum(re,3)./(this.N) );
        end
    end
    
    methods (Access = private)
        function update(this,~,~)
            
            ephi_ = zeros(this.sizeoutC);
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            if numel(this.illu)==1
                illu_ = this.illu.*ones(1, this.Nt);
            else
                illu_ = this.illu;
            end
            
            for nt = 1:this.Nt
                %  frequency grid
                v = 1./(Nx_ *  this.dxy) *( [0:ceil( Nx_/2)-1, -floor( Nx_/2):-1]' -    this.dxy *   sin( this.theta(1,nt))/ this.lambda(nt).*Nx_);
                u = 1./(Ny_ *  this.dxy) * ([0:ceil( Ny_/2)-1, -floor( Ny_/2):-1] -     this.dxy *   sin( this.theta(2,nt))/ this.lambda(nt).*Ny_);
                
                
                [mu,mv] =  meshgrid(  u.^2,  v.^2);
                Mesh = mv + mu;
                
                if (max(Mesh(:))>( this.n0/  this.lambda(nt))^2)
                    mod = illu_(nt).* (Mesh<= ( this.n0/  this.lambda(nt))^2);
                    Mesh(~mod)=0.;
                else
                    mod =  illu_(nt);
                end
                switch  this.type
                    case  this.PUPIL % pupil
                        ephi_(:,:,nt) = mod.*(Mesh<=(this.NA/  this.lambda(nt))^2);
                    case  this.AS % Angular spectrum
                        ephi_(:,:,nt) =  mod.*exp(-2i* pi *   this.z(nt).* sqrt((   this.n0/  this.lambda(nt))^2- Mesh));
                    case  this.FEITFLECK
                        ephi_(:,:,nt)=  mod.*exp(-2i* pi *  this.z(nt).* this.lambda(nt) /  this.n0 * Mesh ./ real(1 + sqrt(1 - ( this.lambda(nt)/ this.n0)^2 *Mesh)));
                    otherwise
                        % separable Fresnel function
                        ephi_(:,:,nt) =mod.* exp(-1i* pi *   this.z(nt).*  this.lambda(nt) /  this.n0 .*Mesh);
                end
            end
            this.ephi = ephi_;
            
            
            
        end
        
    end
end
