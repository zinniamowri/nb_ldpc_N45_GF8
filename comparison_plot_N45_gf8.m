%BER Performance of NG-CNV-SF over GF(8) with 8-QAM, N=45, K=15
%nbldpc_45_2_3_gf8.mat
%algorithm1 (Paper name: )
%EMS
%simulation setup of our own

Eb_No_db=[9.5, 10.5, 11.5, 12.5, 13.5, 14.5];

uncoded=[3.989256e-02, 2.567559e-02, 1.465818e-02, 7.274726e-03, 3.090000e-03, 1.061556e-03];

CNV=[9.513712e-03, 2.653124e-03, 4.823268e-04, 5.377661e-05, 5.777778e-06, 2.222222e-07];

Eb_No_ngdsf=9:1:14;
%T=30, eta=5, w=35, flip_num=2;
NGDSF = [2.336283e-02, 8.714006e-03, 2.515366e-03, 4.659757e-04, 7.652939e-05,  3.407445e-06];

%algorithm1=[1.709049e-02, 5.612815e-03, 1.827029e-03, 4.800646e-04, 1.035000e-04,2.033333e-05];

%data from plotdizitizer
x_for_algo1= [10.5, 11.5, 12.5, 13.38, 14.3, 15.16];
algorithm1=[0.002513489112177349, 0.0009204260882470623, 0.000337055044249586, 0.00003780613901032568, 0.0000056684275986545415, 2.8463849130042736e-7];

%algorithm2=[ 0.01167, 0.00344, 0.000951, 0.0002241, 1.1e-05, 6.333333e-06]

%data from plotdizitizer
x_for_algo2=[3.76, 4.29, 4.58, 5.04, 5.37, 5.8, 6.19];
algorithm2=[0.00005167689353077303, 0.00001809746429760075, 0.0000037089593846081383, 0.0000017754460466292086, 4.6514678944760467e-7, 4.879391338903215e-8, 1.3668915917771727e-8];

x_for_algo3=[6.9, 7.9, 8.8, 9.7, 10.5, 11.5, 12.5];
algorithm3=[0.006564089821917644, 0.0011252435634589665, 0.000575954017022894, 0.0000789781443059708, 0.00001730721332830597, 0.0000021226080457660694, 2.0364055215803732e-7];

%Eb_No_ems=[3.6, 3.8, 4, 4.2, 4.4];
%ber_ems=[3.4327e-05, 2.12528e-05, 7.33333e-06, 4.66667e-06, 1.33333e-06];


figure;

semilogy(Eb_No_db, uncoded, 'rx-', 'LineWidth', 1.2);
grid on;
hold on;

semilogy(Eb_No_db, CNV, 'bx-','LineWidth', 1.2);

semilogy(Eb_No_ngdsf, NGDSF, 'kx-','LineWidth', 1.2);

semilogy(x_for_algo1, algorithm1, '-s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);


semilogy(x_for_algo2, algorithm2, '-s', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);

semilogy(x_for_algo3, algorithm3, '-s','Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2);

%semilogy(Eb_No_ems, ber_ems, 'mx-', 'LineWidth', 1.2);
hold off;

 ylim([10e-9 10e-2]);
 xlim([1 20]);
 xticks([1 5 10 15 20]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 %title('BER comparison for NB-LDPC code (45,2,3) over GF(8)');
 hold off;
 legend( '8-QAM (uncoded)', 'proposed NG-CNV-SF', 'NGDSF','algorithm1','algorithm2','algorithm3', 'Location', 'northeast', 'FontSize', 7);
