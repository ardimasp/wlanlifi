% close all
clear all
clc

f_cut_off = 20e+6;
t = 0:1e-9:200*1e-9;

load('Run1')

f = linspace(0,20e+6,5000);
j = sqrt(-1);
H_LED = 1./(1+j*f/f_cut_off);
H_VLC = zeros(1,length(f));

count = 1;

for fx = f
    
    for t = 1:length(averun2)
        H_VLC(1,count) = H_VLC(1,count) + averun2(t)*exp(-sqrt(-1)*2*pi*fx*(t-1)/1e+9);
    end
    count = count + 1;
    
    
end

H_VLCeff(1,:) = H_VLC(1,:).*H_LED;


plot(f/1e+6,pow2db(abs(H_VLCeff).^2),'linewidth',2)
grid on
xlabel('Frequency [MHz]')
ylabel('|H_{LED}(f)|^{2} [dB]')

hold on
