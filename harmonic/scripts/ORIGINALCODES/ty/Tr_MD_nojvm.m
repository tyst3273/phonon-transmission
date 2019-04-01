@clear
CS='a'; % the prefix of file names 
dT=20 % temperature difference used in NEMD simulations

% Chunk index
kvals=1:20; % devide the whole velocity trajectory into 20 pieces, that is to get 20 independent transmission functions for averages
kn=max(kvals);
tn=50000; % The number of velocity frames used for each piece: tn*kn should <= the total number of velocity frames collected
Nt=tn;

file_forces='Fij.dat'; % the force constant file
file_compact_vels='vels.compact.dat'; % the compactified velocity file

dt_md=0.5e-15; % MD time step, should be changed according to your MD simulation settings
fprintf('Using MD time step %e!\n',dt_md);

kb=1.38063e-23; % boltzmann constant, used for unit conversion
k_tol=1e-5; % Parameter used for creating sparse force constant matrices, useful if the number of atoms is huge.

fid1=fopen(file_forces,'r');

% File format
% NL $NL
% NR $NR
% "atom_id atom_type" of interface atoms, NL+NR lines
% HSTEP $hstep

ss=textscan(fid1,'%s%d',1);
NL=ss{2};
ss=textscan(fid1,'%s%d',1);
NR=ss{2};

ss=textscan(fid1,'%d%d',NL+NR);
atom_ids=ss{1};
atom_types=ss{2};

ids_L=find(atom_types==1);
ids_R=find(atom_types==2);

ss=textscan(fid1,'%s%f',1);
HSTEP=ss{2};

if 1 % Read forces from file
    
    Kii=zeros(3*NL,3*NL);
    Kij=zeros(3*NL,3*NR);

    for i=1:3*NL % Loop over the particles
        if (NL>1e3 && mod(i,1e2)==0)
           fprintf('i=%d/%d\n',i,3*NL); 
        end

        % ITEM: ENTRIES index c_atomids[1] c_atomids[2] c_forces
        ss=textscan(fid1,'%d%f%f%f',NL+NR,'headerlines',10);
        
        % Positive direction shift for the given composite index i
        % (particle number and component)
        inds=ss{1};
        fxs1=ss{2};
        fys1=ss{3};
        fzs1=ss{4};
        
        ss=textscan(fid1,'%d%f%f%f',NL+NR,'headerlines',10);      
        
        % Negative direction shift for the given composite index i
        
        inds=ss{1};
        fxs2=ss{2};
        fys2=ss{3};
        fzs2=ss{4};
        
        % Forces on ids_L and ids_R due to atom displacements of type ids_L
        
        Fii=[fxs1(ids_L)-fxs2(ids_L),fys1(ids_L)-fys2(ids_L),fzs1(ids_L)-fzs2(ids_L)];
        Fij=[fxs1(ids_R)-fxs2(ids_R),fys1(ids_R)-fys2(ids_R),fzs1(ids_R)-fzs2(ids_R)];
        
        Fii=Fii';
        Fij=Fij';
        Fii=Fii(:);
        Fij=Fij(:);

        Kij(i,:)=Fij;
        Kii(i,:)=Fii;
 
    end
end

% The spring constant matrices Kij=d^2V/du_idu_j=-dF_j/du_i come with minus

Kii=-Kii/(2*HSTEP);
Kij=-Kij/(2*HSTEP);

fprintf('Sparsifying with tolerance %d.\n',k_tol);

% Sparse spring constant matrices
Kii_s=Kii;
Kii_s(abs(Kii_s)<k_tol)=0;
Kii_s=sparse(Kii_s);

Kij_s=Kij;
Kij_s(abs(Kij_s)<k_tol)=0;
Kij_s=sparse(Kij_s);

fclose(fid1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NL=length(ids_L);
NR=length(ids_R);

exit_flag=0;
for k=kvals
    if(exit_flag)
        break;
    end
    
    file_vels=strcat(file_compact_vels);
    fprintf('Chunk %d/%d.\n',k,kn);
    
    if (k==1)
        fprintf('Opening the file %s.\n',strcat(file_vels));
        fid2=fopen(file_vels,'r');
    end

    % Atoms NP
    if (k==1)
        ss=textscan(fid2,'%s%d',1);
        NP=ss{2}
        if (NP~=(NL+NR))
		    disp('Mismatch in particle numbers!')
		    return
        end
        % dt_timestep
        ss=textscan(fid2,'%s%d',1);
        d_timesteps=ss{2}
        dt=double(d_timesteps)*dt_md
        % Read the pair ids
        ss=textscan(fid2,'%s%s',1); % "Atom ids:"

        ss=textscan(fid2,'%d',NP);
        ids=ss{1};
        ss=textscan(fid2,'%s',1); %"---------"
    end

 
    disp('Reading velocities...')
    data_v=textscan(fid2,'%.12f',Nt*NP*3);
    data_v=data_v{1};

    if (length(data_v)/(3*NP)~=Nt)
        Nt=length(data_v)/(3*NP);
        fprintf('Changed Nt to %d.\n',Nt);
        if(k>1)
            break; 
        else
            exit_flag=1;
        end
    end

    data_v=reshape(data_v,3*NP,Nt);   

    if (k==kn)
        fclose(fid2);
    end
    vels=data_v;

    vels_fft=fft(vels,[],2)*dt;

    Nf=size(vels_fft,2);

    oms_fft=(0:Nf-1)/(Nf*dt)*2*pi;
    dom=oms_fft(2)-oms_fft(1);

    DoS=abs(vels_fft).^2./repmat(mean(data_v.^2,2),1,Nf)/(Nf*dt);
    Dos_ave = mean(DoS,1); Dos_ave = Dos_ave/sum(Dos_ave);
    
    % Left interface particle indices are in ids_L (from Kij file)
    % Right interface particles are in ids_R (from Kij file)

    vel1_fft=zeros(3*NL,size(vels_fft,2));
    vel2_fft=zeros(3*NR,size(vels_fft,2));

    vel1_fft(1:3:end,:)=vels_fft(3*ids_L-2,:);
    vel1_fft(2:3:end,:)=vels_fft(3*ids_L-1,:);
    vel1_fft(3:3:end,:)=vels_fft(3*ids_L,:);

    vel2_fft(1:3:end,:)=vels_fft(3*ids_R-2,:); 
    vel2_fft(2:3:end,:)=vels_fft(3*ids_R-1,:);
    vel2_fft(3:3:end,:)=vels_fft(3*ids_R,:);

    %Nf_L=size(vel1_fft,2);
    %Nf_R=size(vel2_fft,2);

    data_vL=data_v(1:NL*3,:);
    data_vR=data_v(NL*3+1:end,:);
    
    Dos_L=abs(vel1_fft).^2./repmat(mean(data_vL.^2,2),1,Nf)/(Nf*dt);
    Dos_L_avg = mean(Dos_L,1); Dos_L_avg = Dos_L_avg/sum(Dos_L_avg);

    Dos_R=abs(vel2_fft).^2./repmat(mean(data_vR.^2,2),1,Nf)/(Nf*dt);
    Dos_R_avg = mean(Dos_R,1); Dos_R_avg = Dos_R_avg/sum(Dos_R_avg);

    %% MAY HAVE TO CHANGE mean(data_v,**) to avg across V's in R/L blocks
    %% respectively
  
%     clear vels_fft;

    %%
    Jom=zeros(1,Nf);

    for ki=1:Nf
        if(mod(ki,10000)==0)
            fprintf('k=%d/%d\n',ki,Nf);
        end
        Jom(ki)=-transpose(vel1_fft(:,ki))*Kij_s*conj(vel2_fft(:,ki))/oms_fft(ki);

    end
    
    Jom(1)=0;
    Jom=-sqrt(-1)*(Jom)/(Nf*dt)*2; % Keep the real part
    % Units of K are eV/A^2 and of velocity A/ps (100 m/s) in metal units
    Jom=Jom*1.602e-19/(1e-20)*1e4;
    Jom=Jom/kb/dT;
    if (k==1)
       Jom_raw=zeros(kn+1,Nf);
       Jom_raw(1,:) = oms_fft/2/pi*1e-12;
    end
    Jom_raw(k+1,:) = Jom;
end % Loop over chunks
%%

save(['Tr_', CS, '_raw'],'Jom_raw');

Jom_ave = mean(Jom_raw(2:end,:),1);

%here do gaussian convolution for the raw transmission data
%two convolution width is selected and choose a reasonable one for publication production
dom=oms_fft(2)-oms_fft(1);
dwin=0.3e12*(3*pi); % THz change this value to adjust the convolution width

win=round(dwin/dom);
g = gausswin_my(win); % <-- this value determines the width of the smoothing window
g = g/sum(g);
Tr1=conv(Jom_ave,g,'same');
Dos_g=conv(Dos_ave, g, 'same'); % convolution for Dos using same width as Tr
Dos_L_g=conv(Dos_L_avg, g, 'same');
Dos_R_g=conv(Dos_R_avg, g, 'same');

dwin=0.1e12*(2*pi); % THz
win=round(dwin/dom);
g = gausswin_my(win); % <-- this value determines the width of the smoothing window
g = g/sum(g);      
Tr2=conv(Jom_ave,g,'same');

OM = oms_fft/2/pi*1e-12;
data=[OM; Jom_ave; Tr1; Tr2; Dos_ave];
data = real(data);
fid = fopen(['Tr_', CS],'w');
fprintf(fid,'%10.5f\t%f\t%f\t%f\t%f\n',data);

%figure(1)

%plot(OM, data(3,:), 'k'); %smoothed Tr

%%%%%% Settings for Tr plot %%%%%%%%%%%%%%%%%%%%%%
%xlabel('Frequency (THz)');
%ylabel('Transmission');
%title('Transmission', 'FontSize', 10);
%axis([0 20 0 100]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(2)

%plot(OM, Dos_R_g, 'b', OM, Dos_L_g, 'r'); %Plot local DOS

%%%%%% Settings for DoS plots %%%%%%%%%%%%%%%%%%%%
%xlabel('Frequency (THz)');
%ylabel('Number of States');
%yticks([]);
%title('Density of States', 'FontSize', 10);
%axis([0 20 0 8e-4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Other Useful Plots %%%%%%%%%%%%%%%%%%%%%%%%
%plot(OM, Dos_L_g, 'k'); %smoothed Left Dos_ave
%plot(OM, Dos_L, 'k'); %Left Dos_ave
%plot(OM, Dos_R_g, 'k'); %smoothed Right Dos_ave
%plot(OM, Dos_R, 'k'); %Right Dos_ave
%plot(OM, Dos_g, 'k'); %smoothed Dos_ave
%plot(OM, data(5,:), 'k'); %raw Dos_ave

%plot(OM, data(4,:), 'k'); %raw Tr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear data_v vel1_fft vel2_fft vels data_v

%%%trapz(OM,(-data(3,:)*kb/1.981e-17))*10e12%%%%
%%% unsure of 10e12 factor here. Should be 1e12 right? From from THz to HZ ... 
