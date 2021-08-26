% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper Ha21 on stability of Delay-DAEs

E = [1.0000  0; 0 0]
A = [0.01000  0; 0.60000  -1.0000]
Ad = [-0.2000 -1.0000; -0.0100  -0.8000]

tau = 1.5; h = [0 tau];
tds = tds_create({E},0,{A, Ad},h,'neutral')

tds = tds_create({A, Ad},h)

options=tdsrootsoptions;
v = compute_roots_DDAE(tds,options)

vv = [v.l0; v.l1];
max(real(vv))

figure(1); clf;
subplot(2,2,1)
plot(vv,'r+')
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

print('Example_NavS_Ha21','-depsc')
!epstopdf Example_NavS_Ha21.eps