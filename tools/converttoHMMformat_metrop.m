%function converttoHMMformat()
%
%Converts data saved in masterscriptV2 or V3 to the format Martin uses to
%load HMM code.  The input 'savename' is 'DNA_conc'.
%
%Updated 10/2011 so that this conversion can be done on Metropolis--in
%particular, no uigetdir.  The second input is now a full path to the
%mother directory of data to analyze.
%
%Steph 3/2011

function converttoHMMformat_metrop(savename, dir)

%Have user pick data to convert
%dir = uigetdir('/Volumes/dumbo-4/stephj/TPM data analysis/SJLacI','Choose directory with nolacdata.mat etc');

savedir = fullfile(dir,'HMMformat');
mkdir(savedir)

nolac = load(fullfile(dir,'nolacdata.mat'));
lac = load(fullfile(dir,'lacdata.mat'));
%lacconcat = load(fullfile(dir,'lacdataconcat.mat')); %Don't need this
Gfit = load(fullfile(dir,'Single bead analysis','GaussFit.mat'));
Gfit = Gfit.GaussFit;

totsets = length(nolac.nolacANAL_RMS);

keepbdcount = 1; %Indexes the beads that are going in the final file
Gbdcount = 1; %Indexes the Gaussian fit beads

for i=1:totsets
    
    thisnolacX = nolac.nolacX{i};
    thisnolacY = nolac.nolacY{i};
    thisnolacname = nolac.nolacnames{i};
    thislacX = lac.lacX{i};
    thislacY = lac.lacY{i};
    thislacname = lac.alllacnames{i};
    
    %MasterscriptV3 saves a field called cutpts that records if part or all of
    %a trajectory should be excluded from analysis.  V2 doesn't save this file
    %so the user has to input that information manually:
    if ~isfield(lac,'cutpts')
        thiscutpts = input(strcat('Enter cutpoints for set',int2str(i),'(row vector, keep all = ',int2str(size(thislacX,1)*size(thislacX,2)),'):'));
        cutpts{i} = thiscutpts;
        %MasterscriptV2 and V3 also load the raw data differently so need
        %to reshape if dealing with V2 data:
        tempnolacX = zeros(size(thisnolacX,3),size(thisnolacX,1)*size(thisnolacX,2));
        tempnolacY = zeros(size(thisnolacX,3),size(thisnolacX,1)*size(thisnolacX,2));
        templacX = zeros(size(thislacX,3),size(thislacX,1)*size(thislacX,2));
        templacY = zeros(size(thislacX,3),size(thislacX,1)*size(thislacX,2));

        for r=1:size(thisnolacX,3)
           tempnolacX(r,:) = reshape(transpose(thisnolacX(:,:,r)),size(thisnolacX,1)*size(thisnolacX,2),1);
           tempnolacY(r,:) = reshape(transpose(thisnolacY(:,:,r)),size(thisnolacX,1)*size(thisnolacX,2),1);
           templacX(r,:) = reshape(transpose(thislacX(:,:,r)),size(thislacX,1)*size(thislacX,2),1);
           templacY(r,:) = reshape(transpose(thislacY(:,:,r)),size(thislacX,1)*size(thislacX,2),1);
        end
        clear thisnolacX thisnolacY thislacX thislacY
        thisnolacX = tempnolacX;
        thisnolacY = tempnolacY;
        thislacX = templacX;
        thislacY = templacY;
        clear tempnolacX tempnolacY templacX templacY
        
    else
        thiscutpts = lac.cutpts{i};
    end
    
    if length(lac.lacX) > totsets %Means there was at least one second set
        thislacX2 = lac.lacX{totsets+i};
        thislacY2 = lac.lacY{totsets+i};
        thislacname2 = lac.alllacnames{totsets+i};
        if ~isfield(lac,'cutpts')
            if ~isempty(thislacX2)
                thiscutpts2 = input(strcat('Enter cutpoints for set',int2str(i),'_2 (row vector, keep all = ',int2str(size(thislacX2,1)*size(thislacX2,2)),'):'));
                cutpts{i+totsets} = thiscutpts2;

                templacX2 = zeros(size(thislacX2,3),size(thislacX2,1)*size(thislacX2,2));
                templacY2 = zeros(size(thislacX2,3),size(thislacX2,1)*size(thislacX2,2));

                for r=1:size(thislacX2,3)
                   templacX2(r,:) = reshape(transpose(thislacX2(:,:,r)),size(thislacX2,1)*size(thislacX2,2),1);
                   templacY2(r,:) = reshape(transpose(thislacY2(:,:,r)),size(thislacX2,1)*size(thislacX2,2),1);
                end
                clear thislacX2 thislacY2
                thislacX2 = templacX2;
                thislacY2 = templacY2;
                clear templacX2 templacY2
            end
        else
            thiscutpts2 = lac.cutpts{i+totsets};
        end
    end
    
    for j = 1:size(thislacX,1) %Iterate through the beads in this set
        if thiscutpts(j)~=0 || (exist('thiscutpts2','var') && ~isempty(thiscutpts2) && thiscutpts2(j)~=0) %Otherwise this bead was discarded entirely
            if Gfit(Gbdcount).Approved %Otherwise this bead was discarded during Gauss fitting
                calibration_filename{keepbdcount} = strcat(savename,'_bead',int2str(keepbdcount),'_cal_raw.mat');
                x = [thisnolacX(j,:)',thisnolacY(j,:)'];
                name = strcat(thisnolacname,'_Bead',int2str(j));
                save(fullfile(savedir,calibration_filename{keepbdcount}),'x','name')
                clear x name
                
                x = [thislacX(j,1:thiscutpts(j))', thislacY(j,1:thiscutpts(j))'];
                name = strcat(thislacname,'_Bead',int2str(j));
                if length(lac.lacX) == totsets || isempty(thislacX2) %no second set
                    looping_filename{keepbdcount} = {strcat(savename,'_bead',int2str(keepbdcount),'_area1_trj_raw.mat')};
                    save(fullfile(savedir,looping_filename{keepbdcount}{1}),'x','name');
                    clear x name
                else
                    %It's possible that all data on this bead was discarded
                    %for the first set
                    if thiscutpts(j)~=0
                        trjnames1 = strcat(savename,'_bead',int2str(keepbdcount),'_area1_trj_raw.mat');
                        save(fullfile(savedir,trjnames1),'x','name');
                        clear x name
                    end
                    %It's also possible there was data on this bead during
                    %the first set, but not the second
                    if thiscutpts2(j)~=0
                        trjnames2 = strcat(savename,'_bead',int2str(keepbdcount),'_area2_trj_raw.mat');
                        x = [thislacX2(j,1:thiscutpts2(j))', thislacY2(j,1:thiscutpts2(j))'];
                        name = strcat(thislacname2,'_Bead',int2str(j));
                        save(fullfile(savedir,trjnames2),'x','name');
                        clear x name
                    end
                    if thiscutpts(j)~=0 && thiscutpts2(j)~=0
                        looping_filename{keepbdcount} = {trjnames1,trjnames2};
                    elseif thiscutpts(j)~=0
                        looping_filename{keepbdcount} = {trjnames1};
                    else
                        looping_filename{keepbdcount} = {trjnames2};
                    end
                end
                
                clear trjnames trjnames1 trjnames2
                keepbdcount = keepbdcount+1;
            end
            Gbdcount = Gbdcount+1;
        end
    end
    
    clear thiscutpts thiscutpts2 thislacX thislacX2 thislacY thislacY2 thislacname thislacname2 thisnolacX thisnolacY thisnolacname
    
end

save(fullfile(savedir,'filenames'),'calibration_filename','looping_filename')
if ~isfield(lac,'cutpts')
    save(fullfile(savedir,'cutpts'),'cutpts')
end

