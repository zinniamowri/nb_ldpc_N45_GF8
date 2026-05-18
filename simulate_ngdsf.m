%ngdsf gives good performance for this code and tuning parameters
%code nbldpc_45_2_3_gf8.mat qam8

clear
rng(0)

Eb_No_db =9:1:14;

T=30; %max decoder iteration

eta=5; 
w=35; 

flip_num=2;

p=3; % number of bits per symbol
q = 2^p; 

load('arith_8.mat');
load('nbldpc_45_2_3_gf8.mat'); 

[G, K, pivotCols, freeCols, H_rref] = gf8_generator_from_H(h, add_mat, mul_mat, div_mat);


%h = full(h); % H matrix
N = size(h,2); %length of codeword
M = size(h,1); %number of parity checks
%K=N-M; % msg length
R=K/N; % code rate

[M, ~] = size(h);
str_cn_vn = cell(M,1);

for i = 1:M
    str_cn_vn{i} = find(h(i,:) ~= 0);
end

CN_lst = str_cn_vn;

info_seq=randi([0 q-1], 1, K);

code_seq = gf_vec_mat_mul(info_seq, G, add_mat, mul_mat);

Syndromes = decod_prod(code_seq, h, CN_lst, mul_mat, add_mat);

%QAM-8 mapping
[qam8, qam_binary_map] = generate_qam8_map();


avg_pow = qam8'*qam8/q; %sum of squared magnitudes of all symbols(the total power)/q.
nrm_fct=sqrt(avg_pow); %normalization factor
gf8 = (0:q-1); %GF field symbols
%alph_bin =  logical(fliplr(dec2bin(gf16, p) - 48)); % symbols in binary

%initializing vector size

y = zeros(1, N);
hard_d_cmplx = zeros(1, N);
hard_d_gf16  = zeros(1, N);


FE=zeros(length(Eb_No_db),1);
genFrame=zeros(length(Eb_No_db),1);
iters_cnr=zeros(length(Eb_No_db),1);
BE_Coded=zeros(length(Eb_No_db),1);
BE_unCoded=zeros(length(Eb_No_db),1);


targetFE=500; %maximum FER to be observed
max_gen=1e5; % maximum number of frame to be generated


for i = 1:length(Eb_No_db)


     while(FE(i) < targetFE && genFrame(i)<max_gen)

        genFrame(i)=genFrame(i)+1;

        c(1,1:N) = qam8(code_seq'+1,1); % codeword in complex
        avg_symbol_energy = 1;

        Eb_No_linear(i)= 10.^(Eb_No_db(i)/10);
        No = avg_symbol_energy ./ (p * Eb_No_linear(i)*R); %noise spectral density
        sigma0 = sqrt(No/2)*nrm_fct ; %noise standard deviation
        nse_std=eta*sigma0; %noise perterbation used in flipping function E

        n = sigma0*randn(1,N)+sigma0*randn(1,N)*1i; %complex noise
        y = c + n; % codeword+noise, channel information

        for j=1:N
            distance=abs(qam8-y(j));
            [~,min_idx]= min(distance);
            hard_d_cmplx(j) = qam8(min_idx); % hard decision in complex
            hard_d_gf8(j) = gf8(min_idx);
        end
        
        errors = hard_d_cmplx ~= c; 
        n_errors_hard =sum(errors); %total number of symbol errors after hard decision made
              
        errors_uncoded_bit = zeros(1, K);

        %bit error calculation in hard decision
        for e = 1 : K          
                if hard_d_gf16(e)~=code_seq(e) 
                    s1 = dec2bin(code_seq(e),p);
                    s2 = dec2bin(hard_d_gf8(e),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_uncoded_bit(e) = errors_uncoded_bit(e) + num_diff_bit;                
                end
        end
        
        un_bit_error = sum(errors_uncoded_bit); % no of bit error in each frame
        BE_unCoded(i)=BE_unCoded(i)+un_bit_error; % total no. of bit error in total frame generated in current Eb/No
        
        %calling the decoing function    
        [seqgf,failed,l]= decodeNGDSF(code_seq,hard_d_cmplx, hard_d_gf8,...
        qam8, gf8,y, h, N, M, T, w,add_mat,mul_mat,div_mat,...
        CN_lst, nse_std,qam_binary_map,flip_num,No);  
               
        iters_cnr(i)=iters_cnr(i)+l;

        %bit error calculation for decoded sequence 
    
        errors_coded_bit = zeros(1, K);
            for g = 1 : K            
                if seqgf(g)~=code_seq(g)
                    s1 = dec2bin(code_seq(g),p);
                    s2 = dec2bin(seqgf(g),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_coded_bit(g) = errors_coded_bit(g) + num_diff_bit;                
                end     
            end

        bit_error = sum(errors_coded_bit); % no. of bit error in each frame
        BE_Coded(i)=BE_Coded(i)+bit_error; % total no. of bit error in total frame generated in the current Eb/No

        %iters_cnr(i) = iters_cnr(i)+l;

        if (failed>0)
            FE(i)=FE(i)+1; 
        end       
     end
     fprintf('Eb/No = %.1f dB: BER (coded) = %.6e, BER (uncoded) = %.6e\n', ...
        Eb_No_db(i), ...
        BE_Coded(i) / (genFrame(i) * K * p), ...
        BE_unCoded(i) / (genFrame(i) * K * p));
end

BERunCoded= BE_unCoded ./(genFrame * K *p); %bit error rate for uncoded

BERCoded = BE_Coded ./ (genFrame * K *p); %bit error rate for coded


 figure;
 
 semilogy(Eb_No_db, BERCoded, 'gx-', 'LineWidth', 1.2);
 grid on;
 hold on;

 semilogy(Eb_No_db, BERunCoded, 'rx-','LineWidth', 1.2);
 hold off;

 ylim([10e-8 10e-1]);
 xlim([5 25]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 %title('BER curve LDPC code (45,2,3) over GF(16)');
 hold off;
legend( 'NGDSF', 'Uncoded','Location', 'northeast');











