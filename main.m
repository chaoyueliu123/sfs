
clc;clear;tic


%% par
Nparfor=60;
delay=31;g=0; tao=3;D=0; V_DS=45.36;gg=0.03;
matrix = [1, 0; 0, 0];


delaymin=0; delaymax=119; delaystep=(delaymax-delaymin)/(Nparfor-1);
delay=delaymin:delaystep:delaymax;

%%
for ii = 1:1
    parfor jjj = 1:Nparfor
            [r_fir(ii,jjj,:),n_spike(ii,jjj,:),output_V(ii,jjj,:,:),cISI1(ii,jjj)] = trial(D,V_DS,gg,g,tao,delay(jjj),matrix);  %%
    end
end
%%
for kkk=1:length(cISI1)
    clear aaa bbb ddd
    if kkk==1
        aaa=delay(kkk);
        bbb=cell2mat(cISI1(kkk));
        clear ISI1
        ISI1(2,:)=bbb;
        ISI1(1,:)=aaa;        
    end
    if kkk>1        
        bbb=horzcat(ISI1(2,:),cell2mat(cISI1(kkk)));
        ddd(1:length(cell2mat(cISI1(kkk))))=delay(kkk);
        aaa=horzcat(ISI1(1,:),ddd);
        clear ISI1
        ISI1(2,:)=bbb;
        ISI1(1,:)=aaa;
    end
end

% 创建包含数据的向量

data1 =ISI1(1,:) ;
data2 =ISI1(2,:);
      

% 计算均值和标准差
mean_value = mean(data2);
std_deviation = std(data2);

% 计算CV
CV = (std_deviation / mean_value) * 100;

% 显示结果
fprintf('数据集的均值为 %.2f\n', mean_value);
fprintf('数据集的标准差为 %.2f\n', std_deviation);
fprintf('数据集的CV为 %.2f%%\n', CV);


%% 
       
 for mmm=1:length(cISI1)
    clear AAA BBB DDD
        AAA=delay(mmm);
        BBB=cell2mat(cISI1(mmm));   
        mean_value(mmm) = mean(BBB);
        std_deviation(mmm) = std(BBB);
        CV(mmm) = (std_deviation(mmm) / mean_value(mmm)) * 100;
 end       
        
toc


