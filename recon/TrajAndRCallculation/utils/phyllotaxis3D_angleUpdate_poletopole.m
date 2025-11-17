function [m_adPolarAngle_out m_adAzimuthalAngle_out x y z] = phyllotaxis3D_angleUpdate_poletopole(m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf, nPC)
%%% same function as in IDEA > Phyllotaxistrajectory.h
%%% for angleUpdate functionality just add the uptade of the counter
%%% in idea the trajectory array is (re-)calculated nPC times
%%% here we now also calculate it nPC times and stote it in an nPC times
%%% longer arrray

NProj = m_lNumberOfFrames * m_lProjectionsPerFrame ; % shots * segments (this should not contain nPCs, since the CISSloop is just a repetition of this code. Otherwise trajectory for polar angle chances(
lTotalNumberOfProjections = NProj;   % For UTE, where we are going from pole to pole of the sphere, lTotalNumberOfProjections = NProj

m_adAzimuthalAngle=zeros(1,NProj);
m_adPolarAngle=zeros(1,NProj);

x = zeros (1, NProj);
y = zeros (1, NProj);
z = zeros (1, NProj);

if flagSelf
    N = lTotalNumberOfProjections - m_lNumberOfFrames;
else
    N = lTotalNumberOfProjections ;
end
kost = pi/(2*sqrt(N));

Gn = (1 + sqrt(5))/2;
Gn_ang = 2*pi - (2*pi / Gn);
%Gn_ang = (2*pi / Gn);

for lpc = 1:nPC % outer loop = CISSloop in this case the number of phase cycles nPC

    % count = 1; % original
    %     count = 1 + (m_lProjectionsPerFrame * m_lNumberOfFrames) * (lpc-1); % includes nPC counter
    count = 1 + N * (lpc-1); % correct: includes N instead and nPC counter

    for lk = 1:m_lProjectionsPerFrame % segments
        for lFrame = 1:m_lNumberOfFrames % shots

            linter = lk + (lFrame-1) * m_lProjectionsPerFrame;

            if flagSelf && lk == 1

                m_adPolarAngle(linter) = 0;
                m_adAzimuthalAngle(linter) = 0;

            else


                % pole-to-pole phyllotaxis ESP ' similar update as in IDEA, this counter shouldnt update for the polar angle
                if (count - (N * (lpc-1)) )<=(N/2) % adatp this line in IDEA
                    m_adPolarAngle(linter) = kost * sqrt(2*(count - (N * (lpc-1)) )); % adapt to N in IDEA
                else
                    m_adPolarAngle(linter) =  pi - kost * sqrt(2*(N - (count - (N * (lpc-1)) ))); % adapt N in IDEA
                end

                m_adAzimuthalAngle(linter) = mod ( (count)*Gn_ang, (2*pi) );
                count = count + 1;

            end

            xtmp(linter)= sin(m_adPolarAngle(linter))*cos(m_adAzimuthalAngle(linter));
            ytmp(linter)= sin(m_adPolarAngle(linter))*sin(m_adAzimuthalAngle(linter));
            ztmp(linter)= cos(m_adPolarAngle(linter));


        end
    end
    % count
    % store for all nPCs
    x(:,:,lpc) = xtmp;
    y(:,:,lpc) = ytmp;
    z(:,:,lpc) = ztmp;
    m_adPolarAngle_out(:,:,lpc) = m_adPolarAngle;
    m_adAzimuthalAngle_out(:,:,lpc) = m_adAzimuthalAngle;

end

% close all
% % debug
% figure;
% for ii=1:5
%     subplot(1,5,ii); plot3(x(:,1:100,ii),y(:,1:100,ii),z(:,1:100,ii)); axis off
% end

end