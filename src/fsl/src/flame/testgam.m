a = 0.01;
b = 1/0.01;

cd /usr/local/fsl/src/bpm
dos(sprintf('/usr/local/fsl/src/bpm/testgam %4.4f %4.4f',a,b));

%%% rnd

load /usr/local/fsl/src/bpm/gamrnds -ASCII;

x = [0.1:0.1:20];

figure(2);
clf;
plot(x,gampdf(x,a,1/b)/sum(gampdf(x,a,1/b)))
hold on;
plot(x,gammapdf(a,b,x)/sum(gammapdf(a,b,x)),'y--')

h = hist(gamrnds,x);
plot(x,h/sum(h),'+');

%%% pdf

load /usr/local/fsl/src/bpm/gampdfs -ASCII;

x = [0.1:0.1:20]; 

figure(2);
clf;
plot(x,gampdf(x,a,1/b)/sum(gampdf(x,a,1/b)))
hold on;

plot(x,gampdfs/sum(gampdfs),'+')


