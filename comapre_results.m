%%%%%%
%this script was created in order to plot the  different mesh measure (Y)
% and scenarios measure difference in order to check results consistency.
%%%%%


clear all
close all
% Show differences to homogeneous results
data_homo= toastReadVector('Brest_homo.dat');
data_LESION = toastReadVector('Brest_hetero.dat');
size (data_homo)
% size (data_LESION)
figure;
   plot(real (data_homo));
 hold on
 plot( real(data_LESION));
 xlabel('Detector index: S1 only then S2 only ');
ylabel('Intensity');
%legend('Difference between medium  with random parameter w/without lession',...
legend( ' Intensity measurement for tissue without lession', ' Intensity measurement for tissue with lession');
figure;
 plot((real (data_homo)-real(data_LESION))./real(data_homo));
%  hold on
%  plot(imag (data_homo-data_LESION)./imag(data_homo));

xlabel('Detector index: S1 only then S2 only ');
ylabel('Intensity');
%legend('Difference between medium  with random parameter w/without lession',...
legend( ' real  of the difference',  'imag of the difference','diff' ) %'1 LESION',
%axis equal tight;










% data_homo_680= toastReadVector('breast_dim2_4_homo_680.dat');
% 
% data_close_scat= toastReadVector('breast_dim2_3_les_highscat.dat');
% data_lesion_scatt_rand= toastReadVector('breast_dim2_3_les_highscat_rand.dat');
%  data_1_LESION = toastReadVector('breast_dim2_2_rand.dat');
% % data_LESION = toastReadVector('breast_dim2_3_le_highabs.dat');
%  data_ALL_LESION = toastReadVector('breast_dim2_3_les_rand.dat');
%  data_small_lesion= toastReadVector('breast_dim2_4_les_highscat_homo.dat');
%  data_far__scat_lesion= toastReadVector('breast_dim2_5_les_highscat_homo_850.dat');
%  data_small_lesion_680= toastReadVector('breast_dim2_4_les_highscat_homo_680.dat');
%  data_shape1 = toastReadVector('test.dat');
%  data_shape2 = toastReadVector('test2.dat');


% % size (data_ALL_LESION)
% % data_no_LESION_rand = toastReadVector('breast_all_LESION_rand.dat');
% % data_LESION_rand = toastReadVector('breast_all_LESION_rand_rand.dat');
% 
% %'shape1 no LESION hom','shape1 with LESION ',...
% %     'shape1 with LESION w max para','shape1 with LESION low scattering',...
% %     'shape1 with LESION low absorption','shape1 with LESION w min para'
% % data_LESION_hom = toastReadVector('shape1_LESION_hom.dat');
% % data_no_LESION_1_hom = toastReadVector('shape1_no_LESION_hom.dat');
% % data_no_LESION_wmaxpara = toastReadVector('shape1_LESION_hom_1.dat');
% data_no_LESIONlowscattering = toastReadVector('shape1_LESION_hom_2.dat');
% % data_no_LESION_low absorption = toastReadVector('shape1_LESION_hom_4.dat');
% % data_no_LESIONwminpara = toastReadVector('shape1_LESION_hom_5.dat');
% data_no_LESION_lowscattering_far = toastReadVector('shape1_LESION_hom_far.dat');


%  plot(real(data_1_LESION));
% % hold on
% plot(real(data_ALL_LESION));
% % hold on

% plot(real(data_close_scat));
%  hold on
% %plot((data_close_scat_lesion));
% 
% 
%  plot(real(data_shape1));
%  hold on
%  plot(real(data_shape2));
% 
% xlabel('Detector index S1 only then S2 only ');
% ylabel('Intensity');
% %legend('Difference between medium  with random parameter w/without lession',...
%    legend( ' breast with lesion',  'trapez with lesion',' triang with lesion' ) %'1 LESION',
% %axis equal tight;
% colorbar
% 
% figure
% hold on
% for i=1:size(dlogY,2)
%     ywrap = [dlogY(i:end,i); dlogY(1:i-1,i)];
%     plot(angle,ywrap,'o-');
% end
% axis([0 360 -1 0.1]);
% xlabel('angular source-detector separation');
% ylabel('log intensity perturbation');
