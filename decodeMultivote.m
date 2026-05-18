function [seqgf,failed_init,l]= decodeMultivote(code_seq, hard_d_cmplx, hard_d_gf16,...
    qam16, gf16,y, h, N, M, T, w, add_mat,mul_mat,div_mat,...
    CN_lst, nse_std,qam_binary_map,flip_num, No)

	d_deco_0 = hard_d_cmplx;
    seqgf= hard_d_gf16;
    d_dec_f =hard_d_cmplx;
	l = 0;
    
    d_gf16 = hard_d_gf16;

	S = decod_prod(hard_d_gf16, h, CN_lst, mul_mat, add_mat);
    Sb = double(S==0); % satisfied syndrome S=0 is now considered as 1
	num_satisfied = sum(Sb);
    failed_init= M-num_satisfied;
    % Initialize tracking variables for best codeword
    min_failed = failed_init;
 
    best_d_gf16 = d_gf16;
    
    Hb=double(h>0);

    % Build mapping: GF symbol value 0..15 -> index in gf16/qam16 arrays (1..16)
    gf_to_idx = zeros(1,8);
    for k = 1:8
        gf_to_idx(double(gf16(k)) + 1) = k;
    end

    % VN -> list of neighboring CNs (build once)
    VN_to_CN = cell(1, N);
    for i = 1:M
        vlist = CN_lst{i};   % variables connected to check i
        for t = 1:length(vlist)
            v = vlist(t);
            VN_to_CN{v} = [VN_to_CN{v}, i];
        end
    end

    gf_add = @(a,b) add_mat(a+1, b+1);
    gf_mul = @(a,b) mul_mat(a+1, b+1);
    gf_div = @(a,b) div_mat(a+1, b+1);   

    while (l < T)
         % 1) recompute syndrome for current d_gf16
          
          S = decod_prod(d_gf16, h, CN_lst, mul_mat, add_mat); %new
          Sb = double(S==0); % satisfied syndrome S=0 is now considered as 1
          failed_init = M - sum(Sb);

           % 2) stopping condition
            if failed_init == 0
                best_d_gf16 = d_gf16;
                break;
            end
        
            l = l + 1;

          WSH = w*Sb*Hb;

          %No=2*sigma^2, total noise variance in QAM16
          

          E = (-abs(y - d_dec_f).^2/No) + WSH + nse_std*randn(1,N);


          % Select most unreliable symbols, by given number flip_num
          [~, idx] = mink(E, flip_num);

          temp_gf16 = d_gf16; 

          for jj = 1:length(idx)
                v = idx(jj);
            
                % votes from unsatisfied neighboring checks
                cn_list = VN_to_CN{v}; %neighboring check nodes of variable node v
                votes = [];
            
                %only unsatisfied checks are allowed to vote
                for t = 1:length(cn_list)
                    i = cn_list(t);
            
                    if S(i) == 0
                        continue; % satisfied check gives no vote
                    end

                    %compute what symbol v should be to satisfy check i
            
                    % compute partial_parity_sum = sum_{k in N(i)\v} h(i,k)*x(k)
                    partial_sum = 0; 
                    vlist = CN_lst{i}; %all variable nodes attached to check i
            
                    for u = 1:length(vlist)
                        k = vlist(u);
                        if k == v
                            continue; %skip the symbol we are trying to solve for
                        end
            
                        hik = double(h(i,k));          % get the GF coefficient on edge (i,k)
                        xk  = double(temp_gf16(k));    % get current GF symbol at position k
            
                        term = gf_mul(hik, xk);        %compute hik and sk in GF(16)
                        partial_sum  = gf_add(partial_sum, term);      %accumulate the sum in GF(16)
                    end

                    %Solve for the symbol vote
                    %computes the symbol value that check i wants at position v
            
                    hiv = double(h(i,v));  % coefficient for the edge between check i and variable v
                    vote_sym = gf_div(partial_sum, hiv);   % the symbol value that would make check i satisfied,
                                                   % assuming all the other neighboring symbols stay fixed.
            
                    votes(end+1) = vote_sym;
                end
            
                if isempty(votes)
                    continue;
                end
            
                % pick the most common vote
                % if tie then one closest to y(v) in Euclidean distance is
                % selected
                uniq = unique(votes);
                counts = arrayfun(@(s) sum(votes == s), uniq);
                maxc = max(counts);
                tied_syms = uniq(counts == maxc);
            
                if numel(tied_syms) == 1
                    chosen_sym = tied_syms(1);
                else
                    % tie-break using channel closeness to received y(v)
                    bestd = inf;
                    chosen_sym = tied_syms(1);
                    for q = 1:numel(tied_syms)
                        s = tied_syms(q);
                        qam_idx = gf_to_idx(s + 1);   
                        d = abs(y(v) - qam16(qam_idx))^2;
                        if d < bestd
                            bestd = d;
                            chosen_sym = s;
                        end
                    end
                end

                % assign chosen symbol (convert 0..15 -> gf element via gf16 table)
                temp_gf16(v) = gf16(gf_to_idx(chosen_sym + 1));
            end

            S_temp=decod_prod(temp_gf16, h, CN_lst, mul_mat, add_mat);
            Sb_temp = double(S_temp==0); % satisfied syndrome S=0 is 1
	        num_satisfied_temp = sum(Sb_temp);
            failed_temp = M-num_satisfied_temp;

           if failed_temp < failed_init
                %fprintf("Updating failed_init: %d -> %d\n", failed_init, failed_temp);
                d_gf16 = temp_gf16;             
                failed_init = failed_temp;
                        
                % Save best codeword so far
                if failed_temp < min_failed
                    min_failed = failed_temp;
                    %best_d_gf16 = temp_gf16; %change
                    best_d_gf16 = d_gf16;

                    %fprintf("\nbest codeword updated!\n");    
                end
            end  
    end

    seqgf = best_d_gf16;

    % start of postprocessing     

    % Compute syndrome for the decoded sequence
    S = decod_prod(seqgf, h, CN_lst, mul_mat, add_mat);
    Sb = double(S == 0);  % 1 = satisfied, 0 = unsatisfied

    % Identify unsatisfied check nodes
    unsatisfied_indices = find(Sb == 0);

    % Get variable nodes connected to these unsatisfied check nodes
    VN_candidates = [];
    for i = 1:length(unsatisfied_indices)
        ii = unsatisfied_indices(i);
        VN_list = CN_lst{ii};
        VN_candidates = [VN_candidates, VN_list];
    end
    VN_candidates = unique(VN_candidates);  % remove duplicates

    % Start with best codeword as current decoded codeword
    best_seq = seqgf;
    best_failed = M - sum(Sb);

    for i = 1:length(VN_candidates)
        vn_idx = VN_candidates(i);
        current_symbol = seqgf(vn_idx);
        
        for test_sym = 0:7
            if test_sym ~= current_symbol
                temp_seq = seqgf;
                temp_seq(vn_idx) = test_sym;
                S_check = decod_prod(temp_seq, h, CN_lst, mul_mat, add_mat);
                Sb_check = double(S_check == 0);
                failed_now = M - sum(Sb_check);
    
                if failed_now < best_failed
                    best_failed = failed_now;
                    best_seq = temp_seq;
                end
            end
        end
    end

seqgf = best_seq;
end