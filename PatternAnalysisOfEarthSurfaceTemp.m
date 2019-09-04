datax=VarName1;
datay=VarName2;
%datay_mean=datay;
datay_mean=datay-mean(datay);
class(datax)
dataxx=datetime(datax);
class(dataxx)
%fs=1000;
figure;
plot(dataxx,datay)
datetick('x', 'keepticks','keeplimits')
title('actual plot of temp in years')

Y = fft(datay_mean);
L = length(datay_mean); 

T = 30*24*60*60 
Fs = 1/(T);

f_adj1=(0:L-1)/L;
figure;
plot(f_adj1,abs(Y/L))                          % goes from 0 to 500, with amp of 0 to 1
title('Double Sided Amplitude Spectrum of X(t)')
xlabel('normalized frequency')
ylabel('|P1(f)|')

f_adj2=(0:L-1)*Fs/L;
figure;
plot(f_adj2,abs(Y/L))                          % goes from 0 to 500, with amp of 0 to 1
title('Double Sided Amplitude Spectrum of X(t)')
xlabel('Cyclical frequency')
ylabel('|P1(f)|')

figure;
pwelch(datay_mean,[],[],[],Fs)          %% x which is 0.265 if multiplied by 1000 gives our actual frequency which you can verify with above fft

gmoddd=modwt(datay);
mragggg=modwtmra(gmoddd);

figure;                                      
plot(dataxx(1:300,:),datay_mean(1:300,:), dataxx(1:300,:),mragggg(3,1:300))
datetick('x', 'keepticks','keeplimits')
title('only 300 point comparison of data realized by wavelet')

figure;                                     
plot(dataxx(:,:),datay_mean(:,:), dataxx(:,:),mragggg(3,:))
datetick('x', 'keepticks','keeplimits')
title('Time series Full data realised by the use of wavelet')



figure
findpeaks(mragggg(3,1:500),'MinPeakDistance',10);
datetick('x', 'keepticks','keeplimits')
[ppp,loc]=findpeaks(mragggg(3,:),'MinPeakDistance',2);
%avgtempyearpeak= years(mean(diff(ppp)))                  % so the next temp rise is every 5.1956 years

avgmeantempcycle=years(mean(diff(loc)))

yearnn=year(dataxx);
yearnnn=unique(yearnn);

% for k=1:length(loc)
%     i(k)=loc(k);
% end
% 
% for z=1:length(loc)
%     dataxnew(z)= datax(i(z));
% end    
% dataxnewx=datetime(dataxnew);
% 
% figure
% plot(dataxnewx,ppp);
% datetick('x', 'keepticks','keeplimits');
% title('listing the peaks height and the year at which they occured')

figure;
spectrogram(ppp,flattopwin(266),10,[],Fs);
% spectrogram(ppp,flattopwin(266),10,[],266);
title('spectrogram');

ppp_round=round(ppp,2);
[c,ia,ic]=unique(ppp_round);
ppp_trans=transpose(ppp_round);
%comb= [yearnnn,ppp_trans,ic];
c_tran= transpose(c);
[cc,ian,icc]=unique(ic);

for i=1:length(cc)
    count=0;
    for j=1:length(icc)
        if(cc(i)==icc(j))
            count=count+1;
            appearance_count(i)=count;
        end
    end
end   
appearance_count_tran=transpose(appearance_count);
new_comb=[cc,c_tran,appearance_count_tran];

loc_index=1;
for date_i=1:length(new_comb(:,3))
    if (new_comb(date_i,3) > 1)
         repmore_da(loc_index)=new_comb(date_i,1);
         repmore_val(loc_index)=new_comb(date_i,2);
         repmore_rep_rate(loc_index)=new_comb(date_i,3);
         loc_index=loc_index+1;
    end    
end 
repmore_da_tran=transpose(repmore_da);
repmore_val_tran=transpose(repmore_val);
repmore_rep_rate_tran=transpose(repmore_rep_rate);

repmore_comb=[repmore_da;repmore_val;repmore_rep_rate]
repmore_comb_tran=transpose(repmore_comb);

datay_rounded=round(datay,1)
jjj=1;
for zz=1:length(repmore_da_tran)
    chk_val=repmore_val_tran(zz);
    kk=1;
    while(kk<length(dataxx))
        if (datay_rounded(kk)==chk_val)
            new_matched_date(jjj)= dataxx(kk);
            new_matched_data(jjj)= chk_val;
            jjj=jjj+1;
        end    
        kk=kk+1;
    end    
end

new_new_matched_date=transpose(new_matched_date);
class(new_new_matched_date)
new_new_matched_data=transpose(new_matched_data);
class(new_new_matched_data)       


conv_date=cellstr(new_new_matched_date);
dat_num= datenum(new_new_matched_date); 
extracted_year= year(new_new_matched_date);

%g1_dummy_year_temp_group=[conv_date,new_new_matched_data];
%g1_table=table(g1_year_temp_group);
%g1_sorted_year_temp_g= sort(g1_table)
%plot(g1_year_temp_group(:,1),g1_year_temp_group(:,2))

%scatter(extracted_year,new_new_matched_data)

% new_new_matched_date1=cellstr(new_new_matched_date);
group_data=findgroups(new_new_matched_data)

g1_year_temp_group= [extracted_year,new_new_matched_data,group_data];

% time= datetime(check.date,'InputFormat','mm/dd/yyyy')

check =table;
check.date=new_new_matched_date;
check.data=new_new_matched_data;
check.group_data=group_data;
figure;
gscatter(g1_year_temp_group(:,1),g1_year_temp_group(:,2),g1_year_temp_group(:,3))
title('years vs temperature')
xlabel('years')
ylabel('temperature')


dif_in_year= diff(check.date);

sorted_rows=sortrows(check);

% figure;
% plot(sorted_rows.date, sorted_rows.data)

[sorted_data]=sorted_rows(:,2);
tab_to_cell=table2array(sorted_data);

figure;
boxplot(g1_year_temp_group(:,1),g1_year_temp_group(:,3),'orientation','horizontal')
title('Grouped Temperature set box plot');

for_scat_plt_1= [g1_year_temp_group(1:25,1),g1_year_temp_group(1:25,3)];
for_scat_plt_2= [g1_year_temp_group(26:49,1),g1_year_temp_group(26:49,3)];
for_scat_plt_3= [g1_year_temp_group(50:72,1),g1_year_temp_group(50:72,3)];
for_scat_plt_4= [g1_year_temp_group(73:90,1),g1_year_temp_group(73:90,3)];
for_scat_plt_5= [g1_year_temp_group(91:110,1),g1_year_temp_group(91:110,3)];
for_scat_plt_6= [g1_year_temp_group(111:133,1),g1_year_temp_group(111:133,3)];
for_scat_plt_7= [g1_year_temp_group(134:150,1),g1_year_temp_group(134:150,3)];

intYear1 = diff(for_scat_plt_1(:, 1));
intYear2 = diff(for_scat_plt_2(:, 1));
intYear3 = diff(for_scat_plt_3(:, 1));
intYear4 = diff(for_scat_plt_4(:, 1));
intYear5 = diff(for_scat_plt_5(:, 1));
intYear6 = diff(for_scat_plt_6(:, 1));
intYear7 = diff(for_scat_plt_7(:, 1));

figure;
boxplot(for_scat_plt_1(:,1),for_scat_plt_1(:,2),'orientation','horizontal')
title(['for tmp:5.0 ','Mean: ' num2str(mean(intYear1)) ' years, Standard deviation: ' num2str(std(intYear1)) ' years']);
figure;
boxplot(for_scat_plt_2(:,1),for_scat_plt_2(:,2),'orientation','horizontal')
title(['for tmp:5.1 ','Mean: ' num2str(mean(intYear2)) ' years, Standard deviation: ' num2str(std(intYear2)) ' years']);
figure;
boxplot(for_scat_plt_3(:,1),for_scat_plt_3(:,2),'orientation','horizontal')
title(['for tmp:5.2 ','Mean: ' num2str(mean(intYear3)) ' years, Standard deviation: ' num2str(std(intYear3)) ' years']);
figure;
boxplot(for_scat_plt_4(:,1),for_scat_plt_4(:,2),'orientation','horizontal')
title(['for tmp:5.3 ','Mean: ' num2str(mean(intYear4)) ' years, Standard deviation: ' num2str(std(intYear4)) ' years']);
figure;
boxplot(for_scat_plt_5(:,1),for_scat_plt_5(:,2),'orientation','horizontal')
title(['for tmp:5.4 ','Mean: ' num2str(mean(intYear5)) ' years, Standard deviation: ' num2str(std(intYear5)) ' years']);
figure;
boxplot(for_scat_plt_6(:,1),for_scat_plt_6(:,2),'orientation','horizontal')
title(['for tmp:5.5 ','Mean: ' num2str(mean(intYear6)) ' years, Standard deviation: ' num2str(std(intYear6)) ' years']);
figure;
boxplot(for_scat_plt_7(:,1),for_scat_plt_7(:,2),'orientation','horizontal')
title(['for tmp:5.6 ','Mean: ' num2str(mean(intYear7)) ' years, Standard deviation: ' num2str(std(intYear7)) ' years']);

% figure;
% boxplot(g1_year_temp_group(:,1),g1_year_temp_group(:,3), 'orientation','horizontal')
% figure;
% boxplot(g1_year_temp_group(:,1),g1_year_temp_group(:,3)=='1', 'orientation','horizontal')
%ntYear = diff(Data(:, 1));

% figure
% gscatter(sorted_rows.date,sorted_rows.data,sorted_rows.group_data)
% plot(check.date,check.data,'r*')
% ax=gca;
% set(ax,'yscale','log')



% sorted_tab=sortrows(check,'date')
% sor_k=1;
% for sor_i=1:length(sorted_tab.group_data)
%     sor_chk=sorted_tab.group_data(sor_k);
%     sor_j=1;
%     while(sor_j<length(sorted_tab.group_data))
%         if(sorted_tab.group_data(sor_j)==sor_chk)
%             
%     end    
% end  

% bar(datenum(sorted_tab.date),sorted_tab.data,'FaceColor',[.50 .20 .30] )
% set(hb(1),'FaceColor','b')
% datetick('x','keeplimits','keepticks')

% bar(unique(group_data),unique(new_new_matched_data));
% legend(unique(group_data))

%  g1=[dat_num,new_new_matched_data];
%  g1_group=findgroups(new_new_matched_data)
%  
%  g1_whole_group_tab=[dat_num,new_new_matched_data,g1_group];
%  split_diff=splitapply(@(dat_num){diff(dat_num)},dat_num,g1_group)
%  
%spectrogram(new_new_matched_data)
%dwt(new_new_matched_data,1:100,'sym2','plot')

 %split_diff_result=splitapply(@diff,dat_num,g1_group)
 
% repmore_comb_tran= transpose(repmore_comb);
% k=1;
% z=1;
% for i=1:length(ic)
%     for j=i+1:length(ic) 
%         if ( ppp_trans(i)==ppp_trans(j))
%              
% %          else
% %              ppp_trans_new(i)=ppp_trans(i);
%         end
%     end
%     z=z+1;
% end    
% for i=1:length(ia)
%     k=ia(i);
%     j=1;
%     while j<=length(ppp_trans)
%         if ((j~=k) && (ppp_trans_new(j)~=0))
%             ppp_trans_new(j)=ppp_trans(j);
%         else
%             ppp_trans_new(j)=0;
%         end        
%         j=j+1;        
%     end    
% end
% 
% ppp_tr_wi_year= [yearnnn,ppp_trans_new,ic];
% 
% 
% 
% 
% % j=1;
% % while j <= length(ia)
% %     k=ia(j)
% %     i_indx=1;
% % %     while i_indx<length(ppp_trans)
% % %         if (i_indx ~= k)
% % %             ppp_trans_new(i_indx)=ppp_trans(i_indx);
% % %         else
% % %             ppp_trans_new(i_indx)=0;
% % %         end
% % %         i_indx=i_indx+1;
% % %     end
% %     
% %     j=j+1;
% % end
% % 
% % % 
% % % figure;
% % % plot(1:266,ic,'r.')
% % % title('plotting ic')
% % 
% % % figure;                                 % this shows wrong figure
% % % plot(1:nfft,datayy, 1:nfft,mrag(3,:))
% % % datetick('x', 'keepticks','keeplimits')
% % 
% % % figure;
% % % plot(datexx,datay, datexx,mrag(3,:))
% % % datetick('x', 'keepticks','keeplimits')
% % % 
% % 
% % % figure;
% % % pwelch(datay_mean,[],[],[],fs)
% % % figure
% % % %spectrogram(datay_mean,fs)
% % % spectrogram(datay_mean,flattopwin(1024),10,[],1024);
% % % figure;
% % % cw2=0;
% % % cw1 = cwt(datay_mean,1:10,'sym2','plot'); 
% % % cw2=max(cw1);
% % % figure;
% % % plot(1:1392,cw2)
% % % findpeaks(cw2)
% % 
% % % g=modwt(datay_mean);
% % % 
% % % mrag=modwtmra(g);
% % % figure;
% % % for i =1:7
% % %     subplot(7,1,i)
% % %     plot(mrag(i,:))
% % % end
% % % figure;
% % % plot(mrag(5,:))
% % % % 
% % % figure
% % % subplot(2,1,1)
% % % plot(1:length(datay),datay_mean)
% % % subplot(2,1,2)
% % % plot(1:length(datay),mrag(5,:))
% % % figure
% % % findpeaks(mrag(12,:))
% % 
% % 
% % 
% % 
% % % dataymean=datay-mean(datay);
% % % L= length(datay);
% % % nfft=2^nextpow2(L);
% % % ff= fft(dataymean,nfft);
% % % plot(1:nfft,ff)
