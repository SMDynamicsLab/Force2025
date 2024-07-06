%% Figures from Versaci & Laje 2024
% Load raw data, save preprocessed data


clear;
close all;

% data paths (set current directory to /analysis)
data_dir = '../data/';

load([data_dir 'data_timing_force.mat']);


NTaps_pre=7;

% force profiles
force = [];
for yeari=1:2
    for attentioni=1:2
        for subjecti=1:NSujetos
            for feedbacki=1:2
                for perti=1:2
                    force(yeari,attentioni,subjecti,feedbacki,perti,:,:,:)=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).force(:,prePertZone,:)*1e4;
                end
            end
        end
    end
end



% force peak
f1 = [];
for yeari=1:2
    for attentioni=1:2
        for subjecti=1:NSujetos
            for feedbacki=1:2
                for perti=1:2
                    for triali=1:NTrials
                        for tapi=1:NTaps_pre
                            f1(yeari,attentioni,subjecti,feedbacki,perti,triali,tapi)=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).F1(triali,tapi)*1e4;
                        end
                    end
                end
            end
        end
    end
end


msj = [];
for yeari=1:2
    for attentioni=1:2
        for subjecti=1:NSujetos
            for feedbacki=1:2
                for perti=1:2
                    for triali=1:NTrials
                        for tapi=1:NTaps_pre
                            y=squeeze(force(yeari,attentioni,subjecti,feedbacki,perti,triali,tapi,:));
                            tapDuration_actual=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).tapDuration(triali,tapi);
                            % Mean Squared Jerk (msj)
                            if tapDuration_actual>=176
                                tapDuration_msj=170; %to avoid noise of the discontinuity near the end of the trucated taps
                            else
                                tapDuration_msj=tapDuration_actual;
                            end
                            y1=diff(y); %first derivate
                            y2=diff(y1); %second derivate
                            y3=diff(y2); %third derivate
                            %Mean squared jerk (msj)
                            temp=y3.^2; %ver Hogan 2007 (carpeta kinematics-trayectories)
                            msj(yeari,attentioni,subjecti,feedbacki,perti,triali,tapi)=sum(temp(1:tapDuration_msj))/(2*tapDuration_msj); %complete tap
                        end
                    end
                end
            end
        end
    end
end



% compute variables
asy=[];
for yeari=1:2
    for attentioni=1:2
        for subjecti=1:NSujetos
            for feedbacki=1:2
                f1_(yeari,attentioni,subjecti,feedbacki)=mean(mean([squeeze(f1(yeari,attentioni,subjecti,feedbacki,1,:,:));squeeze(f1(yeari,attentioni,subjecti,feedbacki,2,:,:))],2));
                stdF1(yeari,attentioni,subjecti,feedbacki)=mean(nanstd([squeeze(f1(yeari,attentioni,subjecti,feedbacki,1,:,:));squeeze(f1(yeari,attentioni,subjecti,feedbacki,2,:,:))],0,2));
                msJ(yeari,attentioni,subjecti,feedbacki)=mean(mean([squeeze(msj(yeari,attentioni,subjecti,feedbacki,1,:,:));squeeze(msj(yeari,attentioni,subjecti,feedbacki,2,:,:))],2));
                temp1=([exp(yeari,attentioni).datos(subjecti,feedbacki,1).asy(:,prePertZone);exp(yeari,attentioni).datos(subjecti,feedbacki,2).asy(:,prePertZone)]);
                stdASY(yeari,attentioni,subjecti,feedbacki)=mean(nanstd(temp1,0,2));
                asy(yeari,attentioni,subjecti,feedbacki)=nanmean(nanmean(temp1,2));
                cvF1INTRA(yeari,attentioni,subjecti,feedbacki)=mean([nanstd(exp(yeari,attentioni).datos(subjecti,feedbacki,1).F1(:,prePertZone),0,2)./mean(exp(yeari,attentioni).datos(subjecti,feedbacki,1).F1(:,prePertZone),2);nanstd(exp(yeari,attentioni).datos(subjecti,feedbacki,2).F1(:,prePertZone),0,2)./mean(exp(yeari,attentioni).datos(subjecti,feedbacki,2).F1(:,prePertZone),2)]);
            end
        end
    end
end


% SDA
T1=squeeze(stdASY(1,1,:,1));
T2=squeeze(stdASY(1,1,:,2));
T3=squeeze(stdASY(1,2,:,1));
T4=squeeze(stdASY(1,2,:,2));
T5=squeeze(stdASY(2,1,:,1));
T6=squeeze(stdASY(2,1,:,2));
T7=squeeze(stdASY(2,2,:,1));
T8=squeeze(stdASY(2,2,:,2));
%Pool Year
STDASY=[[T1;T5],[T2;T6],[T3;T7],[T4;T8]];


% CV f1 Intra
T1=squeeze(cvF1INTRA(1,1,:,1));
T2=squeeze(cvF1INTRA(1,1,:,2));
T3=squeeze(cvF1INTRA(1,2,:,1));
T4=squeeze(cvF1INTRA(1,2,:,2));
T5=squeeze(cvF1INTRA(2,1,:,1));
T6=squeeze(cvF1INTRA(2,1,:,2));
T7=squeeze(cvF1INTRA(2,2,:,1));
T8=squeeze(cvF1INTRA(2,2,:,2));
%Pool Year
CVF1INTRA=[[T1;T5],[T2;T6],[T3;T7],[T4;T8]];


%%Mean squared Jerk - msj - complete tap
T1=squeeze(msJ(1,1,:,1));
T2=squeeze(msJ(1,1,:,2));
T3=squeeze(msJ(1,2,:,1));
T4=squeeze(msJ(1,2,:,2));
T5=squeeze(msJ(2,1,:,1));
T6=squeeze(msJ(2,1,:,2));
T7=squeeze(msJ(2,2,:,1));
T8=squeeze(msJ(2,2,:,2));
%Pool Year
MSJ=[[T1;T5],[T2;T6],[T3;T7],[T4;T8]];


% SDF
T1=squeeze(stdF1(1,1,:,1));
T2=squeeze(stdF1(1,1,:,2));
T3=squeeze(stdF1(1,2,:,1));
T4=squeeze(stdF1(1,2,:,2));
T5=squeeze(stdF1(2,1,:,1));
T6=squeeze(stdF1(2,1,:,2));
T7=squeeze(stdF1(2,2,:,1));
T8=squeeze(stdF1(2,2,:,2));
%Pool Year
STDF1=[[T1;T5],[T2;T6],[T3;T7],[T4;T8]];


% F1
T1=squeeze(f1_(1,1,:,1));
T2=squeeze(f1_(1,1,:,2));
T3=squeeze(f1_(1,2,:,1));
T4=squeeze(f1_(1,2,:,2));
T5=squeeze(f1_(2,1,:,1));
T6=squeeze(f1_(2,1,:,2));
T7=squeeze(f1_(2,2,:,1));
T8=squeeze(f1_(2,2,:,2));
%Pool Year
F1=[[T1;T5],[T2;T6],[T3;T7],[T4;T8]];


% Mean force
force_mean=squeeze(nanmean(nanmean(nanmean(force,7),6),5));
T1=squeeze(force_mean(1,1,:,1,:));
T2=squeeze(force_mean(1,1,:,2,:));
T3=squeeze(force_mean(1,2,:,1,:));
T4=squeeze(force_mean(1,2,:,2,:));
T5=squeeze(force_mean(2,1,:,1,:));
T6=squeeze(force_mean(2,1,:,2,:));
T7=squeeze(force_mean(2,2,:,1,:));
T8=squeeze(force_mean(2,2,:,2,:));
FMEAN = [];
FMEAN(:,:,1)=[T1;T5];
FMEAN(:,:,2)=[T2;T6];
FMEAN(:,:,3)=[T3;T7];
FMEAN(:,:,4)=[T4;T8];




save([data_dir 'F1.mat'], 'F1');
save([data_dir 'STDF1.mat'], 'STDF1');
save([data_dir 'CVF1INTRA.mat'], 'CVF1INTRA');
save([data_dir 'STDASY.mat'], 'STDASY');
save([data_dir 'MSJ.mat'], 'MSJ');
save([data_dir 'FORCE.mat'], 'force');
save([data_dir 'FMEAN.mat'], 'FMEAN');



%%

