different = importdata('feedback_B_move_1_180120_Sctime_2000_reinforce_100_4_No_0_1_0.03_0.3_STDP_LiIP_randSc_random_normal_fft_100000.txt')
out_down = find(different(201:300,2)==1);
out_nodif = find(different(201:300,2)==0);
down_sout = sout(out_down,:);
nodif_sout = sout(out_nodif,:);
pos = mean(mean(down_sout(:,1:50)));
nega = mean(mean(down_sout(:,51:100)));
difpos = mean(mean(nodif_sout(:,1:50)));
difnega = mean(mean(nodif_sout(:,51:100)));


