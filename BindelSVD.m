label=1;
matxi = spconvert(load(strcat('cematlab_imagh_34050_6_',num2str(label),'.mtxx')));
matxr = spconvert(load(strcat('cematlab_realh_34050_6_',num2str(label),'.mtxx')));
shiftroffdiag=2;
ndimo = size(matxr,1);
matx = matxr + 1i * matxi + shiftroffdiag*speye(ndimo);

stvx = load(strcat('cestvech_34050_6_',num2str(label),'.stvx'));
stvx = stvx(:,1)+1i*stvx(:,2);


gxx=2.0082;%2.0088;%2.0087
gyy=2.0061;%2.0057
gzz=2.0023;%2.0021
cfact=1e-6*mean([gxx,gyy,gzz])*9.2731e-21/1.05443e-27;

rngeMHz1=300;
rngeMHz2=300;%500;%omarr goes from -rngeMHz to +rngeMHz
NGRD=128;%number of grid points in f1,f2
omarr=2*pi*linspace(-rngeMHz1,rngeMHz2,NGRD);
omarrG=omarr/cfact; %calculate in Gauss units


specarr=zeros(ndimo,NGRD);

parfor i=1:NGRD
    specarr(:,i) = gmres(matx+1i*speye(ndimo)*omarrG(i),stvx,20,1e-7,100);
end
%%
dlmwrite(strcat('uv_34050_6_',num2str(label),'_K.txt'),specarr, 'delimiter', '\t','precision',12);

s=svd(specarr);
figure;
plot(diff(real(stvx'*specarr)));
figure;
semilogy(s,'.-')
