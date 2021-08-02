function [x1,x2,err] = EM_2(p1,p2,data,length_data)%用于两种高斯模型的分类，p1,p2为高斯模型参数向量，p=[k,mu,sigma],data为输入数据列向量,length_data为数据长度
    lnp1 = 0;
    lnp2 = 0;
    r1 = zeros(1,length_data);
    r2 = zeros(1,length_data);%r1,r2为两种分布的后验概率
    for i = 1:length_data%计算对数似然函数
        lnp1 = lnp1+log(p1(1)*normpdf(data(i),p1(2),p1(3))+p2(1)*normpdf(data(i),p2(2),p2(3)))/log(exp(1));
    end
    for i = 1:length_data
        r1(i) = p1(1)*normpdf(data(i),p1(2),p1(3))/(p1(1)*normpdf(data(i),p1(2),p1(3))+p2(1)*normpdf(data(i),p2(2),p2(3)));
        r2(i) = p2(1)*normpdf(data(i),p2(2),p2(3))/(p1(1)*normpdf(data(i),p1(2),p1(3))+p2(1)*normpdf(data(i),p2(2),p2(3)));
    end
    N1 = sum(r1);
    N2 = sum(r2);
    p1(1) = N1/length_data;
    p1(2) = r1*data/N1;                                                                                                                                                                                                             
    p1(3) = sqrt(r1*((data-p1(2)).^2)/N1);
    p2(1) = N2/length_data;
    p2(2) = r2*data/N2;
    p2(3) = sqrt(r2*((data-p2(2)).^2)/N2);
    for i = 1:length_data%计算对数似然函数
        lnp2 = lnp2+log(p1(1)*normpdf(data(i),p1(2),p1(3))+p2(1)*normpdf(data(i),p2(2),p2(3)))/log(exp(1));
    end
    x1 = p1;
    x2 = p2;
    err = abs(lnp2-lnp1);
end