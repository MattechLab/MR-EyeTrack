function [x, y, z,polarAngle,azimuthalAngle] = computePhyllotaxis (N, nseg, nshot, flagSelfNav, flagPlot, flagRosettaTraj)

if nargin < 6
    flagRosettaTraj = false;
end

% flagSelfNav = 1;
if flagRosettaTraj
    [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_FR (nshot, nseg, flagSelfNav);
else
    [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D (nshot, nseg, flagSelfNav);
end%     [polarAngle azimuthalAngle vx vy vz] = phyllotaxis3Dgb (nshot, nseg);
    
%     [polarAngle azimuthalAngle vx vy vz] = archimedean3D (shots, seg);
%     polarAngle = pi - polarAngle;

%r = (-N/2 : N/2-1);
r = (-0.5 : 1/N : 0.5-(1/N));
azimuthal  = repmat(azimuthalAngle,[N 1]);
polar      = repmat(pi/2-polarAngle,[N 1]);
% polar      = repmat(polarAngle,[N 1]);
R          = repmat(r',[1 nshot*nseg]); 

[x, y, z]    = sph2cart(azimuthal,polar,R);

% x = x + N/2;
% y = y + N/2;
% z = z + N/2;

x = reshape(x, [N, nseg, nshot]);
y = reshape(y, [N, nseg, nshot]);
z = reshape(z, [N, nseg, nshot]);

% figure('name','kspace locations')
% plot3(x(:,1:5),y(:,1:5),z(:,1:5),'.'),title('first 5 interleaves')

% figure('name','kspace locations')
% plot3(x,y,z,'.')

%% ... plot 

% % ---- video recording ----
% writerObj = VideoWriter('trajectory.mp4');
% writerObj.FrameRate = 60;
% open(writerObj);

if flagPlot

    figure
    for shot = 1:nshot%201:
        for seg = 1:nseg
            plot3(squeeze(x(:,seg,shot)),squeeze(y(:,seg,shot)),squeeze(z(:,seg,shot)))
            title(['N = ',num2str(N),'  nseg = ',num2str(nseg),'  nshot = ',num2str(nshot)])
            hold on
    %         plot3(squeeze(x(N,seg,shot)),squeeze(y(N,seg,shot)),squeeze(z(N,seg,shot)),'k*')
    %         hold on
    %         plot3(squeeze(x(1,seg,shot)),squeeze(y(1,seg,shot)),squeeze(z(1,seg,shot)),'r*')

                if seg==1
                   hold on
                   plot3(x(end,seg,shot),...
                         y(end,seg,shot),...
                         z(end,seg,shot),'.-k','linewidth',2)
                else
                   plot3([x(end,seg-1,shot),x(end,seg,shot)],...
                         [y(end,seg-1,shot),y(end,seg,shot)],...
                         [z(end,seg-1,shot),z(end,seg,shot)],'.-k','linewidth',2)
                end
            axis([-.5 .5 -.5 .5 -.5 .5])
    %         axis([0 N-1 0 N-1 N/2 N-1])
%             axis([0 N-1 0 N-1 0 N-1])
            pause(.1)
%             % ---- video recording ----
%             frame = getframe;
%             writeVideo(writerObj, frame);
        end
        pause(.3)
    end
    hold off
end
% ---- video recording ----
% close(writerObj);
%%% -----------------------------

end

function [m_adPolarAngle m_adAzimuthalAngle x y z] = phyllotaxis3D (m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf)

NProj = m_lNumberOfFrames * m_lProjectionsPerFrame;
lTotalNumberOfProjections = NProj;   % For UTE, where we are going from pole to pole of the sphere, lTotalNumberOfProjections = NProj

m_adAzimuthalAngle=zeros(1,NProj);
m_adPolarAngle=zeros(1,NProj);

x = zeros (1, NProj);
y = zeros (1, NProj);
z = zeros (1, NProj);

if flagSelf
	kost = pi / ( 2*sqrt(lTotalNumberOfProjections - m_lNumberOfFrames) );
else
    kost = pi / (2*sqrt(lTotalNumberOfProjections));
end    

Gn = (1 + sqrt(5))/2;
Gn_ang = 2*pi - (2*pi / Gn);
%Gn_ang = (2*pi / Gn);
count = 1;

for lk = 1:m_lProjectionsPerFrame
    for lFrame = 1:m_lNumberOfFrames	
        
        linter = lk + (lFrame-1) * m_lProjectionsPerFrame;
        
        if flagSelf && lk == 1
        
            m_adPolarAngle(linter) = 0;
            m_adAzimuthalAngle(linter) = 0;
        
        else
        
            m_adPolarAngle(linter) = kost * sqrt(count);
            m_adAzimuthalAngle(linter) = mod ( (count)*Gn_ang, (2*pi) );
            count = count + 1;
        
        end
        
        x(linter)= sin(m_adPolarAngle(linter))*cos(m_adAzimuthalAngle(linter));
        y(linter)= sin(m_adPolarAngle(linter))*sin(m_adAzimuthalAngle(linter));
        z(linter)= cos(m_adPolarAngle(linter));
        
    end
end

end

function [m_adPolarAngle m_adAzimuthalAngle x y z] = phyllotaxis3D_FR (m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf)

NProj = m_lNumberOfFrames * m_lProjectionsPerFrame;
lTotalNumberOfProjections = NProj;   % For UTE, where we are going from pole to pole of the sphere, lTotalNumberOfProjections = NProj

m_adAzimuthalAngle=zeros(1,NProj);
m_adPolarAngle=zeros(1,NProj);

x = zeros (1, NProj);
y = zeros (1, NProj);
z = zeros (1, NProj);

if (flagSelf)
    kost = pi / ( sqrt(2 * (lTotalNumberOfProjections - m_lNumberOfFrames)));
else
	kost = pi / ( sqrt(2 * lTotalNumberOfProjections));
end

Gn = (1 + sqrt(5))/2;
Gn_ang = 2*pi - (2*pi / Gn);
%Gn_ang = (2*pi / Gn);
count = 1;

for lk = 1:m_lProjectionsPerFrame
    for lFrame = 1:m_lNumberOfFrames	
        
        linter = lk + (lFrame-1) * m_lProjectionsPerFrame;
        
        if flagSelf && lk == 1
        
            m_adPolarAngle(linter) = 0;
            m_adAzimuthalAngle(linter) = 0;
        
        else
        
            m_adAzimuthalAngle(linter) = mod ( (count)*Gn_ang, (2*pi) );
            if count<(lTotalNumberOfProjections/2)
                m_adPolarAngle(linter) = kost * sqrt(count);
            else
    %             m_adPolarAngle(linter) = pi - (kost * sqrt(lTotalNumberOfProjections-count));     
                m_adPolarAngle(linter) = kost * sqrt(lTotalNumberOfProjections-count);     
            end
            count = count + 1;
        end

        x(linter)= sin(m_adPolarAngle(linter))*cos(m_adAzimuthalAngle(linter));
        y(linter)= sin(m_adPolarAngle(linter))*sin(m_adAzimuthalAngle(linter));
        z(linter)= cos(m_adPolarAngle(linter));

    end
end

end