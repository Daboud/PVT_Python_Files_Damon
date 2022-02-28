close all
clear
format bank
tic

col=13; %number of columns in the code. Should be 9 (3D), or 13 (5D)
n=111; %number of points to interpret curved sections. Setting higher makes smoother curves but takes longer
rows=[];
brkpt=1.0; %1 is full plotting. 0.5 is half plotting. 0.1 is 10% plotting
PD=3; %Plot dimensions

shrow=col-3;

fid = fopen('Hammer_CyHi340_CyBa100_Rx1500_Ry162_Zd1_T0.85.txt');  % open the text file


disp('Loading txt file...')

S = textscan(fid,'%s');   % text scan the data
fclose(fid) ;      % close the file
S = S{1} ;
N = cellfun(@(x)str2double(x), S);  % convert the cell array to double 

M = reshape(N',col,numel(N)/col)'; %Convert the vector of values back into the right array format
M(:,6:7)=-M(:,6:7);
if numel(rows)>0; M=M(1:rows,:); end

%------------------ WWIREFRAME BLACK LINES FOR BEAM MOVEMENT---------------

Ms=cumsum(M); %cumulative sum of PVT file columns
shs=M(:,shrow); %this finds where the shutter is actually open
shs=cumsum(shs);
shs(find(M(:,shrow)~=0))=0;
shs(find(shs>0))=1; %ones and zeros showing which rows open the shutter
hits=find(shs==1); %this is the row numbers where the shutter is open
sum(M(:,1))/60

xpo=0;
ypo=0;
zpo=0;
sz=size(M);
if sz(1)>10000
    n=9;
    disp('n auto set to 9 to reduce time')
end

xpt=zeros(1,1);
ypt=zeros(1,1);
zpt=zeros(1,1);
shutl=zeros(1,1);


figure('units','normalized','outerposition',[0 0 1 1])
acx=zeros(1,sz(1));
acy=zeros(1,sz(1));
acz=zeros(1,sz(1));
xyzl=0;

disp('Acceleration and Polynomial Solutions...')

for L=1:sz(1) %for each row
    
    if sz(1)>10000
        if L/20000==floor(L/20000)
            disp(['Accel Calc ',num2str(L/sz(1)*100),' %'])
        end
    end

    t2=M(L,1); %time duration

    v2x=M(L,3); %X speed and distance
    if L==1; v1x=0;
    else v1x=M(L-1,3);
    end
    distx=M(L,2);

    v2y=M(L,5); %Y speed and distance
    if L==1; v1y=0;
    else v1y=M(L-1,5);
    end
    disty=M(L,4);

    v2z=M(L,7); %Z speed and distance
    if L==1; v1z=0;
    else v1z=M(L-1,7);
    end
    distz=M(L,6);

    c1x=(distx-.5*(v2x+v1x)*t2)/(-1/6*t2^3); %analytical solution to a 2nd degree fit for the motion of the stage. Formula is V=c1t^2+c2*t+v1. This is solved for each of X, Y, Z
    c2x=(v2x-c1x*t2^2-v1x)/t2;
    c1y=(disty-.5*(v2y+v1y)*t2)/(-1/6*t2^3);
    c2y=(v2y-c1y*t2^2-v1y)/t2;
    c1z=(distz-.5*(v2z+v1z)*t2)/(-1/6*t2^3);
    c2z=(v2z-c1z*t2^2-v1z)/t2;

    %if the line solution is flat, don't bother interpreting the curve just
    %plot a straight line
    if c1x==0 & c2x==0 & c1y==0 & c2y==0 & c1z==0 & c2z==0
        nn=2;
    else; 
        nn=n; %number of interpretation pts
    end
    
    xyzl=xyzl+nn;

    tt=linspace(0,t2,nn); %time vector
    xp=[];
    yp=[];
    zp=[];
    
    xp=1/3.*c1x.*tt.^3+1/2.*c2x.*tt.^2+v1x.*tt+xpo; %This is the integral of the 2nd order solution above. This track's the stage's movement position throughout the line
    yp=1/3.*c1y.*tt.^3+1/2.*c2y.*tt.^2+v1y.*tt+ypo;
    zp=1/3.*c1z.*tt.^3+1/2.*c2z.*tt.^2+v1z.*tt+zpo;
    
    
    shutl(xyzl-nn+1:xyzl)=shs(L)*ones(1,numel(xp)); %add ones to the shutter's vector if it's open
    xpt(xyzl-nn+1:xyzl)=xp; %long consecutive vector of the stage positions that accumulates throughout the for loop. These will be plotted finally to make the 3D plot
    ypt(xyzl-nn+1:xyzl)=yp;
    zpt(xyzl-nn+1:xyzl)=zp;

    ax=2*c1x*tt+c2x; %Find the acceleration in XYZ for each step
    acx(L)=max(abs(ax));
    ay=2*c1y*tt+c2y;
    acy(L)=max(abs(ay));
    az=2*c1z*tt+c2z;
    acz(L)=max(abs(az));

    xpo=xp(end); %save current position as the for loop continues
    ypo=yp(end);
    zpo=zp(end);
end
disp('Plotting Wireframe...')
if PD==3;
    plot3(xpt,ypt,zpt,'-k','linewidth',0.5) %plot in 3D
elseif PD==2;
    plot(xpt,ypt,'-k','linewidth',0.5) %plot in 2D
end
hold on
%--------------- END WIREFRAME -----------------------
% vx=c1x*tt.^2+c2x*tt+v1x

% ---------------------- START COLOURED LASED SECTIONS -----------------
disp('Plotting Coloured Lased Sections...')
streams=cumsum([0 diff(shutl)]); %binary vector showing where laser is on
starts=[0 diff(shutl)];  %catches the start/end pts of lasing
start=find(starts==1);

cm=jet; %the colormap used to set the colors of the lines
zs=unique(Ms(:,6));
zmin=min(zs);
zmax=max(zs); %finding the min and max z positions so the colormap uses full range

count=0;

tot=numel(start);
totr10=round(tot,-1)/10;

q=rand(1,numel(start));
p=find(q>brkpt);
start(p)=[];

for i=start
    
    count=count+1;
    if sz(1)>10000
        if (count/totr10)==floor(count/totr10)
            disp(['Colour Plotting ',num2str(count/totr10*10/brkpt),' %'])
        end
    end
    
    nowkill=min(find(starts(i:end)==-1))+i-1;
%     disp(i)
    
    z=zpt(i); %setting colormap based on starting z-position of each lasing segment
    zu=(z-zmin)/(zmax-zmin); %setting the colours of the lines
    cp=round(zu*255)+1;
    if isnan(cp); 
        co=cm(end,:);
    else
        co=cm(cp,:)/1;
    end
    
%     co=cm(end,:);
    
    
    if PD==3;
        plot3(xpt(i:nowkill),ypt(i:nowkill),zpt(i:nowkill),'-','LineWidth',5,'color',co) %plot coloured lines in 3D
    elseif PD==2;
        plot(xpt(i:nowkill),ypt(i:nowkill),'-','LineWidth',5,'color',co) %plot coloured lines in 2D
    end
end


% Plot the accelleration of the XYZ stages throughout on double y-axis plot
figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(acx,'-','linewidth',1.5,'color','g')
hold on
plot(acy,'-b','linewidth',1.5)
ylim([0 max([acx acy])*1.2])
xlabel('step')
ylabel('Acceleration')

if sum(acz~=0)
yyaxis right
plot(acz,'linewidth',1.5)
ylim([0 max([acz])*1.2])
end

legend('X','Y','Z')

figure(1)
xlabel('X')
ylabel('Y')
zlabel('Z')
toc
% 
% Msa=Ms(:,6)-Ms(:,4);
% Msa=round(Msa,4)
% figure
% plot(Msa)
% title('Relative Z Position')
% msau=unique(Msa);









