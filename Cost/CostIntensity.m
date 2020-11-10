classdef CostIntensity < Cost
    % COMCI Library: Intensity cost function
    %
    % Implement the proximity operator for intensity measurement  described in [1]
    % $$ \phi(x) = \sum_n f(\sum_k  |x_{k,n}|^2 - y_n , w) $$ where f is the log-likelihood of 
    % either Normal or Poisson distributions.
    %
    % :param sz: size of the input
    % :param y:  the measured intensity   
    % :param w:  if type=Gaussian, w is the inverse variance of the measurement 
    %            if type=Poisson, b is the background flux.
    % :param index: index along wich dimension the intensity is summed
    % :param type: noise model 'Gaussian' or 'Normal' for Gaussian noise or 'Poisson' for Poisson noise.
    %
    % **References**
    %
    % [1] Ferréol Soulez, Eric Thiébaut, Antony Schutz, André Ferrari, Frédéric Courbin, and Michael Unser, "Proximity operators for phase retrieval," Appl. Opt. 55, 7412-7421 (2016) 
    %
    % **Example** C=...
    %
    % Obj = CostIntensity(sz,y,w,  'Gaussian')
    %
    %
    % See also :class:`Map` :class:`Cost`
    %% Properties
    %
    %%
    
    %     Copyright (C) 2015-2018 F. Soulez ferreol.soulez@epfl.ch
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
    
    properties (SetAccess = protected,GetAccess = public)
        I % Intensity value
        w % Precision (inverse variance) in the Gaussian case (default 1) or  background the Poisson noise case (default 0)
        index % index along wich dimension are computed the finite differences
        kerdims % ker dimensions
        imdims % im dimensions
        szx;    % size of the input
        widx=[]; % index of invalid data
        noisemodel  % flag for noise model:
        %*  0 stationnary Gaussian (default)
        %*  1 non stationnary Gaussian with precision vector w
        %*  2 Poisson with background w=0;
        %*  3 Poisson with background w!=0;
        %*  4 Laplace
        Pcost  %base cost ( log(y!) ) for Poisson
    end
    
    methods
        function this = CostIntensity(sz,y,w,varargin)
            this.name='CostIntensity';

            assert( isnumeric(y),  'y should be numeric');
            
           this.sizein = sz;
           this.isConvex=false; 
            this.y = y(:) ;
            index=[];
            this.noisemodel = 0;
            if nargin >2
                if isscalar(w)
                    this.w = w;
                else
                    if isnumeric(w)
                        assert( isequal(size(w),size(y)),  'w does not have the right size: [%d, %d, %d,%d]',size(y));
                        %   assert(any(w(:)>=0),'weight or background must be positive');
                        
                        this.w = w(:);
                    else
                        error('w should be numeric');
                    end
                end
                
            end
            for c=1:length(varargin)
                if (isvector(varargin{c})&& isnumeric(varargin{c}))
                    index = varargin{c};
                else
                    switch varargin{c}
                        case{'Gaussian','Normal'}
                            if isscalar(this.w)
                                this.noisemodel = 0;
                            else
                                this.noisemodel = 1;
                                this.widx = (this.w>0);
                                this.w = this.w(this.widx);
                                this.y = this.y(this.widx);
                            end
                        case('Poisson')
                            if any(y(:)<0)
                                warning('Intensity should be positive');
                            end
                            this.y = max(y(:),0);
                            
                            if (isscalar(this.w))
                                if(this.w==0)
                                    this.noisemodel = 2;
                                else
                                    this.noisemodel = 3;
                                end
                            else
                                this.widx = (this.w>=0);
                                this.w = this.w(this.widx);
                                this.y = this.y(this.widx);
                                if (all(this.w==0))
                                    this.noisemodel=2;
                                else
                                    this.noisemodel = 3;
                                end
                            end
                            % Approx log(y!)
                            xx = max(this.y-1,1);
                            this.Pcost = sum((xx-1./2.).*log(xx) - xx + 1/2.*log(2.*pi) + 1./(12.*xx));
                            
                        case('Laplace')
                            this.noisemodel = 4;
                    end
                end
            end
            if (~isempty(index))
                assert(isvector(index) ,'The index a vector');
                this.index = sort(index,'descend');
            else
                this.index =[];
            end
            
            
            
            
        end
        
        function cost = applySafe(this,x) % get the function cost
            
            
            Ix = abs(x).^2;
            if (~isempty(this.index))
                for n=this.index
                    Ix = sum(Ix,n);
                end
            end
            if (~isempty(this.widx))
                Ix = Ix(this.widx);
            else
                Ix = Ix(:);
            end
            
            switch(this.noisemodel)
                case {0,1} % Gaussian noise
                    
                    
                    res = this.w.* ( Ix - this.y).^2;
                    cost = sum(res(:));
                case {2,3} % Poisson noise
                    Ix = Ix(:) + this.w(:);
                    t = (this.y==0);
                    cost = this.Pcost;
                    if any(t)
                        cost = cost + sum(Ix(t) );
                    end
                    if any(~t)
                        cost = cost + sum(Ix(~t) -  this.y(~t) .* log(Ix(~t)));
                    end
                case 4 % Laplace noise
                    res = abs( abs(x(:)).^2 - this.y);
                    cost = sum(res(:));
            end
        end
    
        function z = applyProxSafe(this,x,alpha) % Apply the operator
            %   assert(isscalar(alpha),'alpha must be a scalar');
            assert( isnumeric(x),  'x should be numeric');
            
            this.szx = size(x);
            
            ndms = length(this.szx);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.szx(2) ==1
                ndms = 1;
            end
            T = true(ndms,1);
            T(this.index)=false;
            this.kerdims = this.szx;
            this.kerdims(T)=1;
            this.imdims = this.szx;
            this.imdims(~T)=1;
            
            
            z=x;
            x=abs(x);
            %     y = y./x;
            if (~isempty(this.index))
                nx = x.^2;
                for n=this.index
                    nx = (sum(nx,n));
                end
                nx = sqrt(nx);
                x = nx;
            else
                nx= x;
            end
            
            if (~isempty(this.widx))
                ymodt = zeros_(numel(x),1);
                ymodt(~this.widx) =x(~this.widx);
                x = x(this.widx);
                if (~isscalar(alpha))
                    alpha = alpha(this.widx);
                end
            end
            x = x(:);
            ymod = zeros_(numel(x),1);
            
            switch this.noisemodel
                case {0,1} % Gaussian case
                    % Amounts to find the hight real root of
                    % x^3 - 3q x = 2r;
                    q = -1./3* ( 1./(4*alpha * this.w) - this.y);
                    r = -0.5*x./(alpha*4*this.w);
                    p = q.^3 - r.^2;
                    
                    t = (p>=0);
                    % Case p>0 % three real roots
                    if any(t)
                        qt = sqrt(q(t));
                        theta = acos(r(t)./(qt.^3))./3.0;
                        ymod(t)= -2* qt.* min( min(cos(theta),cos(theta+2*pi/3)),cos(theta+4*pi./3));
                        
                    end
                    
                    if any(~t)
                        sr = sign(r(~t));
                        sr(sr==0)=1; % assume 0 is positive
                        pt =  sr.*nthroot( sqrt(-p(~t)) + abs(r(~t)) , 3);
                        ymod(~t) = -pt-q(~t)./pt;
                    end
                    %
                case 2 % Poisson noise without background (solve a quadratic equation)
                    b = alpha*x;
                    ymod = 0.5./(2+ alpha)*(b+sqrt(b.^2+ 8 * (2+ alpha) * this.y));
                    
                case 3 % Poisson noise + background
                    % the solution is the highest real roots of
                    % x^3 + a2 x^2 + a1 x + a0
                    
                    a2 = -1./(1 + 2.*alpha).*x;
                    a1 = this.w - 2.*alpha./(1+2.*alpha) .* this.y;
                    a0 = a2.* this.w;
                    
                    q = -1./9. .*( 3.*a1 - a2.^2);
                    r = -1./54. .* (9.*a1.*a2 - 27.* a0 - 2.* a2.^3);
                    
                    p = q.^3 - r.^2;
                    
                    t = (p>=0);
                    % Case p>0 % three real roots
                    if any(t)
                        qt = sqrt(q(t));
                        theta = acos(r(t)./(qt.^3))./3.0;
                        ymod(t)= -2* qt.* min( min(cos(theta),cos(theta+2*pi/3)),cos(theta+4*pi./3))- 1./3..*a2(t);
                        
                    end
                    
                    if any(~t)
                        sr = sign(r(~t));
                        sr(sr==0)=1; % 0 is positive
                        pt =  sr.*nthroot( sqrt(-p(~t)) + abs(r(~t)) , 3);
                        ymod(~t) = -pt-q(~t)./pt- 1./3..*a2(~t);
                    end
                case 4 % Laplace
                    ymod = sqrt(this.y);
                    st  = sign(x.^2 - this.y);
                    t = (st==1);
                    ymod(t) = max( 1./(2.*alpha+1).*x(t), sqrt(this.y(t)));
                    
                    if alpha<0.5
                        t = (st==-1);
                        ymod(t) = min( 1./(2.*alpha-1).*x(t), sqrt(this.y(t)));
                    end
                    
                    
            end
            
            if (~isempty(this.widx))
                ymodt(this.widx) = ymod;
                ymod= ymodt;
            end
            
            
            znx = (nx~=0);
            if (any(znx))
                ymod(znx)= ymod(znx)./nx(znx);
            end
            z = z.* reshape(repmat(reshape(ymod,this.imdims),this.kerdims),this.szx); % y = x ./ |x| * ymod
            
        end
    end 
    methods (Access = protected)
        function z = applyProx_(this,x,alpha) % Apply the operator
            %   assert(isscalar(alpha),'alpha must be a scalar');
            assert( isnumeric(x),  'x should be numeric');
            
            this.szx = size(x);
            
            ndms = length(this.szx);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.szx(2) ==1
                ndms = 1;
            end
            T = true(ndms,1);
            T(this.index)=false;
            this.kerdims = this.szx;
            this.kerdims(T)=1;
            this.imdims = this.szx;
            this.imdims(~T)=1;
            
            
            z=x;
            x=abs(x);
            %     y = y./x;
            if (~isempty(this.index))
                nx = x.^2;
                for n=this.index
                    nx = (sum(nx,n));
                end
                nx = sqrt(nx);
                x = nx;
            else
                nx= x;
            end
            
            if (~isempty(this.widx))
                ymodt = zeros_(numel(x),1);
                ymodt(~this.widx) =x(~this.widx);
                x = x(this.widx);
                if (~isscalar(alpha))
                    alpha = alpha(this.widx);
                end
            end
            x = x(:);
            ymod = zeros_(numel(x),1);
            
            switch this.noisemodel
                case {0,1} % Gaussian case
                    % Amounts to find the hight real root of
                    % x^3 - 3q x = 2r;
                    q = -1./3* ( 1./(4*alpha * this.w) - this.y);
                    r = -0.5*x./(alpha*4*this.w);
                    p = q.^3 - r.^2;
                    
                    t = (p>=0);
                    % Case p>0 % three real roots
                    if any(t)
                        qt = sqrt(q(t));
                        theta = acos(r(t)./(qt.^3))./3.0;
                        ymod(t)= -2* qt.* min( min(cos(theta),cos(theta+2*pi/3)),cos(theta+4*pi./3));
                        
                    end
                    
                    if any(~t)
                        sr = sign(r(~t));
                        sr(sr==0)=1; % assume 0 is positive
                        pt =  sr.*nthroot( sqrt(-p(~t)) + abs(r(~t)) , 3);
                        ymod(~t) = -pt-q(~t)./pt;
                    end
                    %
                case 2 % Poisson noise without background (solve a quadratic equation)
                    b = alpha*x;
                    ymod = 0.5./(2+ alpha)*(b+sqrt(b.^2+ 8 * (2+ alpha) * this.y));
                    
                case 3 % Poisson noise + background
                    % the solution is the highest real roots of
                    % x^3 + a2 x^2 + a1 x + a0
                    
                    a2 = -1./(1 + 2.*alpha).*x;
                    a1 = this.w - 2.*alpha./(1+2.*alpha) .* this.y;
                    a0 = a2.* this.w;
                    
                    q = -1./9. .*( 3.*a1 - a2.^2);
                    r = -1./54. .* (9.*a1.*a2 - 27.* a0 - 2.* a2.^3);
                    
                    p = q.^3 - r.^2;
                    
                    t = (p>=0);
                    % Case p>0 % three real roots
                    if any(t)
                        qt = sqrt(q(t));
                        theta = acos(r(t)./(qt.^3))./3.0;
                        ymod(t)= -2* qt.* min( min(cos(theta),cos(theta+2*pi/3)),cos(theta+4*pi./3))- 1./3..*a2(t);
                        
                    end
                    
                    if any(~t)
                        sr = sign(r(~t));
                        sr(sr==0)=1; % 0 is positive
                        pt =  sr.*nthroot( sqrt(-p(~t)) + abs(r(~t)) , 3);
                        ymod(~t) = -pt-q(~t)./pt- 1./3..*a2(~t);
                    end
                case 4 % Laplace
                    ymod = sqrt(this.y);
                    st  = sign(x.^2 - this.y);
                    t = (st==1);
                    ymod(t) = max( 1./(2.*alpha+1).*x(t), sqrt(this.y(t)));
                    
                    if alpha<0.5
                        t = (st==-1);
                        ymod(t) = min( 1./(2.*alpha-1).*x(t), sqrt(this.y(t)));
                    end
                    
                    
            end
            
            if (~isempty(this.widx))
                ymodt(this.widx) = ymod;
                ymod= ymodt;
            end
            
            
            znx = (nx~=0);
            if (any(znx(:)))
                ymod(znx)= ymod(znx)./nx(znx);
            end
            z = z.* reshape(repmat(reshape(ymod,this.imdims),this.kerdims),this.szx); % y = x ./ |x| * ymod
            
        end
        function cost = apply_(this,x) % get the function cost
            
            
            Ix = abs(x).^2;
            if (~isempty(this.index))
                for n=this.index
                    Ix = sum(Ix,n);
                end
            end
            if (~isempty(this.widx))
                Ix = Ix(this.widx);
            else
                Ix = Ix(:);
            end
            
            switch(this.noisemodel)
                case {0,1} % Gaussian noise
                    
                    
                    res = this.w.* ( Ix - this.y).^2;
                    cost = sum(res(:));
                case {2,3} % Poisson noise
                    Ix = Ix(:) + this.w(:);
                    t = (this.y==0);
                    cost = this.Pcost;
                    if any(t)
                        cost = cost + sum(Ix(t) );
                    end
                    if any(~t)
                        cost = cost + sum(Ix(~t) -  this.y(~t) .* log(Ix(~t)));
                    end
                case 4 % Laplace noise
                    res = abs( abs(x(:)).^2 - this.y);
                    cost = sum(res(:));
            end
        end
    end
    
end
