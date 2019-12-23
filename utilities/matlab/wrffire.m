classdef wrffire
    
    % WRFFIRE is a collection of functions for plotting and analysing
    % output from WRF-Fire simulations.
    %
    % Trent Mattner
    % November 2012

    properties (Constant)

        gravity = 9.81;

    end
    
    properties
        
        ncid; % file handle
        dx;   % east-west grid spacing
        dy;   % north-south grid spacing
        dt;   % time-step
        nx;   % number of cells in east-west direction
        ny;   % number of cells in north-south direction
        nz;   % number of cells in vertical direction
        nt;   % number of time steps
        nsgx;
        nsgy;
        
    end
    
    methods
        
        function obj = wrffire(filename)
            
            % WRFFIRE creates a handle to the data file and reads the
            % grid parameters.
            
            obj.ncid = netcdf.open(filename,'NC_NOWRITE');
            obj.dx = obj.getAtt('DX');
            obj.dy = obj.getAtt('DY');
            obj.dt = obj.getAtt('DT');
            obj.nx = obj.getDim('west_east');
            obj.ny = obj.getDim('south_north');
            obj.nz = obj.getDim('bottom_top');
            obj.nt = obj.getDim('Time');
            obj.nsgx = obj.getDim('west_east_subgrid');
            obj.nsgy = obj.getDim('south_north_subgrid');
            
        end
        
        function delete(obj)
            
            % DELETE closes the file.
            
            disp('Deleting...')
            netcdf.close(obj.ncid);
            
        end
        
        function attr = getAtt(obj,attrname)
            attr = netcdf.getAtt(obj.ncid,netcdf.getConstant('NC_GLOBAL'),attrname,'double');
        end
        function dim = getDim(obj,dimname)
            try
                id = netcdf.inqDimID(obj.ncid,dimname);
                [~,dim] = netcdf.inqDim(obj.ncid,id);
            catch
                dim = [];
            end
        end
        function var = getVar(obj,varname,itime)
            % Reads in ARW variable named varname. itime is the time index.
            if strfind(varname,'altitude')
                ph = obj.getVar('PH',itime);    % base geopotential
                phb = obj.getVar('PHB',itime);  % perturbation geopotential
                var = (ph + phb)/9.8;           % altitude
            elseif strfind(varname,'pressure')
                var = obj.getVar('P',itime) + obj.getVar('PB',itime);
            elseif strfind(varname,'temperature')
                p = obj.getVar('pressure',itime);
                t = obj.getVar('T',itime);
                p00 = obj.getVar('P00',itime);
                var = (t + 300).*(p/p00).^(2/7) - 273.16;
            elseif strfind(varname,'volume/dewpoint')
                qv = obj.getVar('QVAPOR',itime);
                p = obj.getVar('pressure',itime)/100; % convert to hPa
                tdc = qv.*p./(0.622+qv);       % partial pressure water vapour (hPa)
                tdc = max(tdc,0.001);
                var = (243.5*log(tdc) - 440.8)./(19.48 - log(tdc));
            elseif strfind(varname,'surface/dewpoint')
                qv = obj.getVar('Q2',itime);
                p = obj.getVar('PSFC',itime)/100; % convert to hPa
                tdc = qv.*p./(0.622+qv);       % partial pressure water vapour (hPa)
                tdc = max(tdc,0.001);
                var = (243.5*log(tdc) - 440.8)./(19.48 - log(tdc));
            elseif strfind(varname,'surface/relhumidity') % google, percentage
                qv = obj.getVar('Q2',itime);
                p = obj.getVar('PSFC',itime)/100; % hPa
                tc = obj.getVar('T2',itime) - 273.16;
                e = qv.*p./(0.622+qv);
                es = 6.1121*(1.0007+3.46e-6*p).*exp(17.502*tc./(240.97+tc)); % saturation vapour pressure (Buck, 1996)
                var = 100*e./es;
            elseif strfind(varname,'volume/relhumidity') % google, percentage
                qv = obj.getVar('QVAPOR',itime);
                p = obj.getVar('pressure',itime)/100; % hPa
                tc = obj.getVar('temperature',itime);
                e = qv.*p./(0.622+qv);
                es = 6.1121*(1.0007+3.46e-6*p).*exp(17.502*tc./(240.97+tc)); % saturation vapour pressure (Buck, 1996)
                var = 100*e./es;
            else
                varid = netcdf.inqVarID(obj.ncid,varname);
                [~,~,dimids,~] = netcdf.inqVar(obj.ncid,varid);
                ndim = length(dimids);
                if dimids(ndim) == 0                % variable is time-dependent
                   varsize = zeros(1,ndim);
                   for i=1:ndim
                       [~,varsize(i)] = netcdf.inqDim(obj.ncid,dimids(i));
                   end
                   start = [zeros(1,ndim-1) itime-1];
                   count = [varsize(1:ndim-1) 1];
                   var = netcdf.getVar(obj.ncid,varid,start,count,'double');
                else
                   var = netcdf.getVar(obj.ncid,varid,'double');
                end
            end
        end
        function var = getTimes(obj)
            varid = netcdf.inqVarID(obj.ncid,'Times');
            var = netcdf.getVar(obj.ncid,varid)'; 
        end
    end

    methods(Static)

        function fv = isosurface(x,y,z,u,varargin)

            % ISOSURFACE returns the face, vertex and (optionally) colour data
            % of an isosurface of u defined on the structured but nonuniform 
            % grid given by the 3d coordinate arrays x, y and z. The isosurface 
            % vertices are first found in index space, then they are mapped to 
            % physical space by interpolation. Plot using PATCH.

            fv = isosurface(u,varargin{:});
            if size(fv.vertices,1) ~= 0
                vx = interp3(x,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3));
                vy = interp3(y,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3));
                vz = interp3(z,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3));
                fv.vertices = [vx vy vz];
            end

        end

        function fv = isoclip(fv,xmin,xmax,ymin,ymax,zmin,zmax)

            % ISOCLIP sets vertices outside the box defined by xmin, xmax,
            % ymin, ymax, zmin, zmax to NaNs. PATCH will not plot these.

            if ~isempty(fv.vertices)
                i = fv.vertices(:,1) < xmin | ...
                    fv.vertices(:,1) > xmax | ...
                    fv.vertices(:,2) < ymin | ...
                    fv.vertices(:,2) > ymax | ...
                    fv.vertices(:,3) < zmin | ...
                    fv.vertices(:,3) > zmax;
                j = find(i);
                if ~isempty(j)
                    fv.vertices(j,:) = NaN(length(j),size(fv.vertices,2));
                end
            end

        end   
     
        function [xmin,xmax,ymin,ymax] = bbox(x,y)

            % BBOX returns lateral domain boundaries, defined as the largest
            % box that will fit inside the grid.

            xmin = max(x(1  ,:));
            xmax = min(x(end,:));
            ymin = max(y(:,1  ));
            ymax = min(y(:,end));

        end

        function [x3,y3] = grid3(n,x,y)

            % GRID3 stacks 2d lateral coordinate arrays vertically.

            x3 = repmat(x,[1 1 n]);
            y3 = repmat(y,[1 1 n]);

        end
        
        function [i,j] = coord2sub(x,y,xref,yref)
            
            % COORD2SUB finds the subscripts of the elements in x and y 
            % that are closest to xref and yref.
            
            dist = (x - xref).^2 + (y - yref).^2;
            [~,k] = min(dist(:));
            [i,j] = ind2sub(size(x),k);
            
        end

        function u = interp3d(x3,y3,z3,u3,x,y,z,varargin)
            
            % INTERP3D interpolates the gridded data u3(x3,y3,z3)
            % at the points (x,y,z). It is assumed that x3 and y3
            % are on a regular grid and independent of z3.

            n = length(x);
            nz = size(z3,3);
            zxy = zeros(n,nz);
            uxy = zeros(n,nz);
            for k = 1:nz
                zxy(:,k) = interpn(x3(:,:,1),y3(:,:,1),z3(:,:,k),x(:),y(:),varargin{:});
                uxy(:,k) = interpn(x3(:,:,1),y3(:,:,1),u3(:,:,k),x(:),y(:),varargin{:});
            end
            u = zeros(size(x));
            for i = 1:n
                u(i) = interpn(zxy(i,:),uxy(i,:),z(i),varargin{:});
            end

        end
        
        function ux = interpx(u)
            
            % INTERPX(U) interpolates U to the midpoint in the x direction,
            % assuming constant grid spacing.
            
            ux = 0.5*(u(1:end-1,:,:) + u(2:end,:,:));
            
        end
        
        function uy = interpy(u)
            
            % INTERPY(U) interpolates U to the midpoint in the y direction,
            % assuming constant grid spacing.
            
            uy = 0.5*(u(:,1:end-1,:) + u(:,2:end,:));
            
        end
        
        function un = interpn(u,varargin)
            
            % INTERPN(U) interpolates U vertically from cell faces to 
            % cell centres on a nonuniform vertical grid.  
            
            % INTERPN(U,RDN,RDN_W) interpolates U vertically from
            % cell centres to cell faces on a nonuniform vertical grid.
            % RDN and RDN_W are arrays containing the reciprocal of the
            % vertical grid spacing between cell centres and cell faces,
            % respectively. Note that the output is located on the
            % interior cell-faces only. See equation (3.28), WRF Technical
            % Document.

            if numel(varargin) == 2
                rdn = varargin{1};
                rdn_w = varargin{2};
                un = 0.5*rdn.*(u(:,:,2:end  ).*rdn_w(:,:,1:end-1) ...
                             + u(:,:,1:end-1).*rdn_w(:,:,2:end) );
            else
                un = 0.5*(u(:,:,2:end) + u(:,:,1:end-1));
            end
            
        end
        
        function uz = interpz(u)
            
            uz = 0.5*(u(:,:,1:end-1) + u(:,:,2:end));
            
        end
        
        function dudx = ddx(u,dx)
            
            % DDX(U,DX) calculates the derivative of U in the x direction
            % assuming constant grid spacing DX.
            
            dudx = (u(2:end,:,:) - u(1:end-1,:,:))/dx;
            
        end

        function dudy = ddy(u,dy)
            
            % DDY(U,DY) calculates the derivative of U in the y direction
            % assuming constant grid spacing DY.
            
            dudy = (u(:,2:end,:) - u(:,1:end-1,:))/dy;
            
        end

        function dudn = ddn(u,rdn)
            
            % DDN(U,RDN) calculates the derivative of U in the surface normal
            % direction. RDN is the reciprocal of the wall-normal grid
            % spacing.
            
            dudn = (u(:,:,2:end) - u(:,:,1:end-1)).*rdn;
            
        end
        
        function u3 = scatterxy(u,n)
            
            % SCATTERXY(U,N) takes a vector U and makes N(1)*N(2) copies in 
            % the xy plane.
            
            u3 = reshape(u,1,[]);
            u3 = u3(ones(n(1)*n(2),1),:);
            u3 = reshape(u3,[n(1) n(2) length(u)]);
            
        end
        
        function u3 = scattern(u,n)
            
            % SCATTERN takes a 2d array U and makes N copies in the surface
            % normal direction.
            
            u3 = repmat(u,[1 1 n]);
            
        end
        
        function [wx,wy,wz] = vorticity(u,v,w,gz,dx,dy,rdn,rdn_w,mx,my,mx_u,my_v)
            
            % VORTICITY(...) calculates the vorticity.
            
            nx = size(w,1); % number of cells in x direction
            ny = size(w,2); % number of cells in y direction
            nz = size(u,3); % number of cells in z direction
            
            % Scatter scaling factors vertically
            
            mx_u = wrffire.scattern(mx_u,nz); % (nx+1,ny,nz)
            my_v = wrffire.scattern(my_v,nz); % (nx,ny+1,nz)
            
            % Scale velocities
            
            U = u./mx_u; % (nx+1,ny,nz)
            V = v./my_v; % (nx,ny+1,nz)
            
            % Calculate velocity derivatives at cell vertices
            
            dvdx = wrffire.ddx(V,dx); % (nx-1,ny+1,nz)
            dudy = wrffire.ddy(U,dy); % (nx+1,ny-1,nz)
            dwdx = wrffire.ddx(w,dx); % (nx-1,ny,nz+1)
            dwdy = wrffire.ddy(w,dy); % (nx,ny-1,nz+1)
            
            % Interpolate to cell centres
            
            dvdx = wrffire.interpx(wrffire.interpy(dvdx)); % (nx-2,ny,nz)
            dudy = wrffire.interpy(wrffire.interpx(dudy)); % (nx,ny-2,nz)
            dwdx = wrffire.interpx(wrffire.interpn(dwdx)); % (nx-2,ny,nz)
            dwdy = wrffire.interpy(wrffire.interpn(dwdy)); % (nx,ny-2,nz)
            
            % Remove edge cells
            
            dvdx = dvdx(:,2:end-1,2:end-1); % (nx-2,ny-2,nz-2)
            dudy = dudy(2:end-1,:,2:end-1); % (nx-2,ny-2,nz-2)
            dwdx = dwdx(:,2:end-1,2:end-1);
            dwdy = dwdy(2:end-1,:,2:end-1);
            
            % Scatter vertical grid spacing laterally
            
            rdn_w = wrffire.scatterxy(rdn_w,[nx ny]);    % (nx,ny,nz)
            rdn = wrffire.scatterxy(rdn(2:end),[nx ny]); % (nx,ny,nz-1)
            
            % Gradients of geopotential
            
            dgzdx = wrffire.ddx(gz,dx); % (nx-1,ny,nz+1)
            dgzdy = wrffire.ddy(gz,dy); % (nx,ny-1,nz+1)
            dgzdn = wrffire.ddn(gz,rdn_w); % (nx,ny,nz)
            
            dgzdx = wrffire.interpn(wrffire.interpx(dgzdx)); % (nx-2,ny,nz)
            dgzdy = wrffire.interpn(wrffire.interpy(dgzdy)); % (nx,ny-2,nz)
            
            % Velocity gradients in vertical direction
            
            dudn = wrffire.ddn(wrffire.interpx(U),rdn); % (nx,ny,nz-1)
            dvdn = wrffire.ddn(wrffire.interpy(V),rdn); % (nx,ny,nz-1)
            dwdn = wrffire.ddn(w,rdn_w); % (nx,ny,nz)
            
            % Interpolate to cell centres
                     
            dudn  = wrffire.interpn(dudn); % (nx,ny,nz-2)
            dvdn  = wrffire.interpn(dvdn); % (nx,ny,nz-2)
            
            % Remove edge cells
            
            dgzdx = dgzdx(:,2:end-1,2:end-1);
            dgzdy = dgzdy(2:end-1,:,2:end-1);
            dgzdn = dgzdn(2:end-1,2:end-1,2:end-1);
            dudn = dudn(2:end-1,2:end-1,:);
            dvdn = dvdn(2:end-1,2:end-1,:);
            dwdn = dwdn(2:end-1,2:end-1,2:end-1);
            
            dudz = wrffire.gravity./dgzdn.*dudn;
            dvdz = wrffire.gravity./dgzdn.*dvdn;
            
            mx = wrffire.scattern(mx(2:end-1,2:end-1),nz-2);
            my = wrffire.scattern(my(2:end-1,2:end-1),nz-2);
            
            % Vorticity at last
            
  %          wx =  my.*(dwdy + dgzdy./dgzdn.*dwdn - dvdz);
  %          wy = -mx.*(dwdx + dgzdx./dgzdn.*dwdn - dudz);    
  %          wz =  mx.*my.*(dvdx + dgzdy./dgzdn.*dvdn ...
  %                      - dudy - dgzdx./dgzdn.*dudn);
                     
            wx =  my.*(dwdy - dgzdy./dgzdn.*dwdn - dvdz);
            wy = -mx.*(dwdx - dgzdx./dgzdn.*dwdn - dudz); 
            wz =  mx.*my.*(dvdx - dgzdx./dgzdn.*dvdn ...
                         - dudy + dgzdy./dgzdn.*dudn);
            
        end
        
        
        function skewTlogP()
            % Creates a Skew-T Log-P plot
            
            % Constants
            axisColor = [0.6 0.6 0.3];
            isoColor = [0 0.8 0.8];
            
            % Isotherms
            [t,p] = meshgrid(-100:10:160,logspace(2,log10(1050),20));
            [x,y] = wrffire.skewmap(t,p);
            plot(x,y,'Color',axisColor)
            hold on
            axis equal
            
            % Isobars
            p = 200:100:1000;
            [~,y] = wrffire.skewmap(p,p);
            [x,y] = meshgrid(linspace(-40,0,5),y);
            plot(x',y','Color',axisColor)
            
            % Dry adiabats
            [t,p] = meshgrid(250:10:440,logspace(2,log10(1050),20));
            t = t - 273.16;
            tda = wrffire.tda(t,p);
            [x,y] = wrffire.skewmap(tda,p);
            plot(x,y,'--','Color',isoColor)
            % Labels
            t = 320:10:440;
            p = 105*ones(1,13);
            t = t - 273.16;
            tda = wrffire.tda(t,p);
            [x,y] = wrffire.skewmap(tda,p);
            text(x,y,'320|330|340|350|360|370|380|390|400|410|420|430|440','Color',isoColor)
            
            % Saturated adiabats
            [t,p] = meshgrid(8:4:32,logspace(log10(200),log10(1050),20));
            tsa = wrffire.tsa(t,p);
            [x,y] = wrffire.skewmap(tsa,p);
            plot(x,y,'--','Color',isoColor)
            % Labels
            t = [20 24 28 32];
            p = 190*ones(1,4);
            tsa = wrffire.tsa(t,p);
            [x,y] = wrffire.skewmap(tsa,p);
            text(x,y,'20|24|28|32','Color',isoColor)           
            
            % Mixing ratio
            [w,p] = meshgrid([0.4 1 2 3 4 5 8 12 20],logspace(log10(700),log10(1050),20));
            tmr = wrffire.tmr(w,p);
            [x,y] = wrffire.skewmap(tmr,p);
            plot(x,y,'--','Color',isoColor)
            % Labels
            w = [0.4 1 2 3 4 5 8 12 20];
            p = 680*ones(1,9);
            tmr = wrffire.tmr(w,p);
            [x,y] = wrffire.skewmap(tmr,p);
            text(x,y,'0.4|1|2|3|4|5|8|12|20','Color',isoColor)
            hold off
            
            % Axes
            set(gca,'XColor',axisColor)
            set(gca,'YColor',axisColor)
            pticks = (1000:-100:100);
            [~, yticks] = wrffire.skewmap(zeros(size(pticks)),pticks);
            set(gca,'YTick',yticks);
            set(gca,'YTickLabel','1000|900|800|700|600|500|400|300|200|100')
            paxis = [1050 100];
            [~, yaxis] = wrffire.skewmap(zeros(size(paxis)),paxis);
            tticks = (-40:10:40);
            [xticks, ~] = wrffire.skewmap(tticks,1050*ones(size(tticks)));
            set(gca,'XTick',xticks);
            set(gca,'XTickLabel','-40|-30|-20|-10|0|10|20|30|40')
            taxis = [-35 50];
            [xaxis, ~] = wrffire.skewmap(taxis,1000*ones(size(taxis)));
            axis([xaxis yaxis])
            
            % Block out upper right hand corner
            p = [680 400 90 90 680];
            t = [40 5.07 -43.4231 -25.7434 40];
            [x,y] = wrffire.skewmap(t,p)
            patch(x,y,'w','EdgeColor',axisColor)
            
        end
        function [x,y] = skewmap(t,p)
            % This is the mapping given by Stipanuk
            x = 0.1408*t - 10.53975*log10(p); % + 31.61923;
            y = -11.5*log10(p); % + 34.5;
        end
        function t = tda(ts,p)
            t = (ts+273.16).*(p/1000).^0.288 - 273.16;
        end
        function e = esat(t)
            e = 10.^(23.832241 - 5.02808*log10(t) ...
                - 1.3816e-7*10.^(11.344-0.0303998*t) ...
                + 8.1328e-3*10.^(3.49149 - 1302.8844./t) ...
                - 2949.076./t);
        end
        function t = tsa(ts,p)
            b = -2.6518986;
            t = ts + 273.16;
            w = 622*wrffire.esat(t)./(1000 - wrffire.esat(t));
            a = t./exp(b*w./t);
            t = 253.16;
            d = 120;
            for i = 1:12
                d = d/2;
                w = 622*wrffire.esat(t)./(p - wrffire.esat(t));
                x = a.*exp(b*w./t) - t.*(1000./p).^0.288;
                t = t + sign(x)*d;
            end
            t = t - 273.16;

%             t = 273.16 + ts;
%             w = @(T,P) 622*wrffire.esat(T)./(P - wrffire.esat(T));
%             a = t./exp(b*w(t,1000)./t);
%             f = @(T) a.*exp(b*w(T,p)./T) - T.*(1000./p).^0.288;
%             for i=1:20
%                 fk = f(t);
%                 dfk = (f(t+0.01) - fk)/0.01;
%                 t = t - fk./dfk;
%             end
%             t = t - 273.16;
        end
        function t = tmr(w,p)
            a = 0.0498646455;
            b = 2.4082965;
            c = 280.23475;
            d = 38.9114;
            f = 0.0915;
            g = -1.2035;
            m = w.*p./(622 + w);
            t = 10.^(a*log10(m) + b) - c + d*(m.^f + g).^2;
        end
    end
end
