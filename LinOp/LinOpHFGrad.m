classdef LinOpHFGrad <  LinOp
    %% LinOpHFGrad :
    %  Matlab Linear Operator Library
    %
    
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp DFT Sfft iSFFT
    
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
    
    properties (SetAccess = protected,GetAccess = public)
        mtfx
        mtfy
        ndms
        N
    end
    methods
        function this =  LinOpHFGrad(sz)
            
            this.name ='Linop HFGrad';
            %   this.isinvertible=false;
            %   this.iscomplex = true;
            
            
            this.sizein =sz;
            
            this.sizeout = [sz 2];
            
            this.ndms = length(this.sizein);
            
            
            this.N = prod(sz);
            
            
            der = zeros(sz);
            der(1,1) = 1;
            der(end,1) = -1;
            this.mtfx = fft2(der);
            der = zeros(sz);
            der(1,1) = 1;
            der(1,end) = -1;
            
            this.mtfy = fft2(der);
        end
    end
    methods  (Access = protected)
        function y = apply_(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d]',this.sizein);
            
            y =  cat(3,ifft2( this.mtfx.*x),ifft2( this.mtfy.*x));
            
        end
        function y = applyAdjoint_(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            
            y =  1./(this.N) .*(conj(this.mtfx).*fft2( x(:,:,1)) + conj(this.mtfy).*fft2( x(:,:,2)));
            
        end
        function y = applyHtH_(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y =  1./(this.N) .*((real(this.mtfx).^2 + imag(this.mtfx).^2 + real(this.mtfy).^2 + imag(this.mtfy).^2) .* x);
            
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein,1./   this.N .*(real(this.mtfx).^2 + imag(this.mtfx).^2 + real(this.mtfy).^2 + imag(this.mtfy).^2) );
        end
        
    end
end
