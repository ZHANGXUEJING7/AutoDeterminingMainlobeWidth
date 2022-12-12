

clear;
clc;
close all;

f_size=12;
lambda=1;
ele_num=60;
sensor_pos=[0:(ele_num-1)]'*(lambda/2);

s_angle=-25;
diff_angle=0.2;
seek_theta=[-90:diff_angle:90];
seek_num=length(seek_theta);
center_index=round((s_angle+90)/diff_angle)+1;
st_all=exp(1j*2*pi*sensor_pos*sin(seek_theta/180*pi)/lambda);

t_st=st_all(:,center_index);
resp_qui=abs(t_st'*st_all);
pattern_qui0=20*log10(resp_qui./max(resp_qui)+1.0e-6);

main_left_angle=s_angle-1;
main_left_index=round((main_left_angle+90)/diff_angle)+1;
main_right_angle=s_angle+1;
main_right_index=round((main_right_angle+90)/diff_angle)+1;

slack_width=15;

main_left2_angle=s_angle-slack_width;
main_left2_index_org=round((main_left2_angle+90)/diff_angle)+1;
main_right2_angle=s_angle+slack_width;
main_right2_index_org=round((main_right2_angle+90)/diff_angle)+1;


slack_region=[[main_left2_index_org:main_left_index] [main_right_index:main_right2_index_org]];



pattern_desired=zeros(1,seek_num);

for reg_desired=1:seek_num
    cur_angle=seek_theta(reg_desired);
    if(cur_angle>main_left_angle)&&(cur_angle<main_right_angle)
        pattern_desired(reg_desired)=0/0;
    else
            pattern_desired(reg_desired)=-40;
    end
end
pattern_desired(main_left_index+1:main_right_index-1)=0/0;


opt_diff=1;
slack_num=length(slack_region);


gain_mat=zeros(slack_num-2,slack_num);

penal_vec=zeros(1,slack_num);
for reg=1:(slack_num-2)
    if(reg<=slack_num/2-1)
        gain_mat(reg,reg)=1;
        gain_mat(reg,reg+1)=-1;
        penal_vec(reg)=slack_num-reg;
    else
        gain_mat(reg,reg+1)=-1;
        gain_mat(reg,reg+2)=1;
        penal_vec(reg)=reg;
    end
end

%*********************************
penal_vec(slack_num-1)=slack_num-1;
penal_vec(slack_num)=slack_num;
%*********************************
% penal_vec=ones(1,slack_num);

penal_vec=rand(1,slack_num)*2;

vec_e=zeros(slack_num,1);
vec_e(slack_num/2)=1;
vec_e(slack_num/2+1)=1;
vec_e=vec_e';


main_left2_index=main_left2_index_org;
main_right2_index=main_right2_index_org;


side_lobe_v=10.^(pattern_desired./10);
iteration_index=0;
iteration_MAX=4;


%*********************************
s_threshold=1e-4;

s_absolute_value=zeros(iteration_MAX,seek_num)/0;
penal_absolute_value=zeros(iteration_MAX,seek_num)/0;
penal_normalized_value=zeros(iteration_MAX,seek_num)/0;
%*********************************


while(iteration_index<iteration_MAX)
    iteration_index=iteration_index+1
    cvx_begin
    variable w_sca(ele_num,1) complex
    variable s(slack_num,1)
    minimize penal_vec*s
    
    subject to
    w_sca'*t_st==1;
    s>=0;
    gain_mat*s<=0;
    slack_index=1;
    for loop1=1:opt_diff:seek_num
        if(loop1<=main_left_index)||(loop1>=main_right_index)
            a_cur=st_all(:,loop1);
            if(any(slack_region==loop1))
                abs(w_sca'*a_cur)<=sqrt(side_lobe_v(loop1))+s(slack_index);
                slack_index=slack_index+1;
            else
                abs(w_sca'*a_cur)<=sqrt(side_lobe_v(loop1));
            end
        end
    end
    cvx_end
    
    
    %*********************************
    s_absolute_value(iteration_index,slack_region)=s;
    penal_absolute_value(iteration_index,slack_region)=penal_vec;
    
    penal_vec_max=max(penal_vec);
    penal_vec_normal=penal_vec./penal_vec_max(1);
    penal_normalized_value(iteration_index,slack_region)=penal_vec_normal;
    %*********************************


    %*********************************
    [find_s_small_x,find_s_small_y]=find(s<=s_threshold);
    main_left2_index_new=slack_region(1);
    main_right2_index_new=slack_region(end);
    for s_x_reg=1:length(find_s_small_x)
        if(find_s_small_x(s_x_reg)+slack_region(1)-1<main_left_index)
            main_left2_index_new=find_s_small_x(s_x_reg)+slack_region(1);
        elseif(find_s_small_x(s_x_reg)+slack_region(1)-1>main_left_index)
            main_right2_index_new=slack_region(end)-(length(slack_region)-find_s_small_x(s_x_reg))-1;
            break;
        else
            
        end
    end    
    
    slack_region=[[main_left2_index_new:main_left_index] [main_right_index:main_right2_index_new]];
    
    slack_num=length(slack_region);
    
    s_new=s;
    s_new(find_s_small_x)=[];
    penal_vec=1./(s_new.'+eps);
    %*********************************
    
    
    
    
    %*********************************
    gain_mat=zeros(slack_num-2,slack_num);    
    for reg=1:(slack_num-2)
        if(reg<=slack_num/2-1)
            gain_mat(reg,reg)=1;
            gain_mat(reg,reg+1)=-1;
        else
            gain_mat(reg,reg+1)=-1;
            gain_mat(reg,reg+2)=1;
        end
    end
    %*********************************
    
    s
    s_new
    
    resp_qui=abs(w_sca'*st_all);
    pattern_cvx=20*log10(resp_qui./max(resp_qui)+1.0e-6);
    figure;
    plot(seek_theta,pattern_qui0,'color','c','LineStyle','-.','LineWidth',1.5);hold on;
    plot(seek_theta,pattern_cvx,'LineWidth',2,'color','r','LineStyle','-');hold on;
    plot(seek_theta,pattern_desired,'k--','LineWidth',1);hold off;
    ylim([-70 10]);
    xlim([-90,90]);
    set(gca,'fontsize',f_size);   
    
end


pattern_cvx3=pattern_cvx;







