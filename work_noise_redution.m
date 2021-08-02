L = 256;%÷°≥§
R = 128;%÷°“∆
win = hann(L);%¥∞∫Ø ˝
N = length(data);
k = mod(N,R);
data1 = [data;zeros(L-k,1)];
out = zeros(length(data1),1);
e = 1e-10;%±‹√‚∑÷ƒ∏¡„÷µ
alpha = 0.1;
ng = 0.05;
r = 0.998;
b = 0.8;
H = zeros(L,1);
p = zeros(L,2);
p_min = zeros(L,1);
power_y = zeros(L,1);
snr_prio = zeros(L,1);
snr_post = zeros(L,1);
sr = zeros(L,1);
flag = zeros(L,1);
ap = 0.2;
pos = zeros(L,1);
ad = 0.85;
as = zeros(L,1);
noise_power = zeros(L,1);
for i = 0:(length(data1)/R-2)
    frame_i = data1((i*R+1):(i*R+L)).*win;
    frame_k = fft(frame_i,L);
    power_y = abs(frame_k).*abs(frame_k);
    if(i==0)
        p_min = power_y;
        noise_power = power_y;
        p(:,1) = power_y;
    end
    cos_y = real(frame_k)./abs(frame_k);
    sin_y = imag(frame_k)./abs(frame_k);
    p(:,2) = ng*p(:,1)+(1-ng)*power_y;
    p(:,1) = p(:,2);
    for j = 1:L
        if(p_min(j)<p(j,2))
            p_min(j) = r*p_min(j)+(1-r)/(1-b)*(p(j,2)-b*p(j,1));
        else
            p_min(j) = p(j,2);
        end
    end
    sr = p(:,2)./p_min;
    for k = 1:(3*L/16+1)
        if(sr(k)>2)
            flag(k) = 1;
        else
            flag(k) = 0;
        end
    end
    for k = (3*L/16+2):(L/2+1)
        if(sr(k)>5)
            flag(k) = 1;
        else
            flag(k) = 0;
        end
    end
    flag((L/2+2):L) = flag((L/2):-1:2);
    pos = ap*pos+(1-ap)*flag;
    noise_power = pos.*noise_power+(1-pos).*power_y;
    snr_post = power_y./noise_power;
    snr_prio = alpha*snr_prio+(1-alpha)*((snr_post-1)+abs(snr_post-1))/2;
    H = snr_prio./(1+snr_prio+e);
    power_x = H.*power_y;
    XK = sqrt(power_x).*cos_y+1i*(sqrt(power_x).*sin_y);
    out_data = real(ifft(XK,L));
    out((i*R+1):(i*R+L)) = out((i*R+1):(i*R+L))+out_data;
end
audiowrite('pross_audio1.wav',out,16000);