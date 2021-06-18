function [V,W]=raydec1station(vert, north, east, time, fmin, fmax, fsteps, cycles, dfpar, nwind)

% RAYDEC1STATION(VERT, NORTH, EAST, TIME, FMIN, FMAX, FSTEPS, CYCLES, DFPAR, NWIND) 
%    calculates the ellipticity of Rayleigh waves for the
%    input data VERT, NORTH, EAST and TIME for a single station 
%    for FSTEPS frequencies (on a logarithmic scale) between 
%    FMIN and FMAX, using CYCLES periods for the stacked signal 
%    and DFPAR as the relative bandwidth for the filtering.
%    The signal is cut into NWIN different time windows and 
%    RayDec is applied to each of them. 
% 
%    VERT, NORTH, EAST and TIME have to be data matrices 
%    (N x 1 or 1 x N) of equal sizes
%
%    suggested values: CYCLES = 10
%                      DFPAR = 0.1
%                      NWIND such that the single time windows are about 10 minutes long
%
%    Code written by Manuel Hobiger, 
%    Laboratoire de Géophysique Interne et Tectonophysique (LGIT), Grenoble, France, 2008-10
%    Last change: 18/06/2021

v1=vert;n1=north;e1=east;t1=time;
if size(v1,2)>size(v1,1)           % formatting the signals such that the vectors are N x 1
   v1=transpose(v1);
   n1=transpose(n1);
   e1=transpose(e1);
   t1=transpose(t1);
end

% setting up
K0=size(v1,1);
K=K0/nwind;
tau=t1(2)-t1(1);
DTmax=30;
fnyq=1/(2*tau);
fstart=max(fmin,1/DTmax);
fend=min(fmax,fnyq);
flist=zeros(fsteps,1);
constlog=(fend/fstart)^(1/(fsteps-1));
fl=fstart*constlog.^(cumsum(ones(fsteps,nwind))-1);
el=zeros(fsteps,nwind);

% loop over the time windows
for ind1=1:nwind
    vert=detrend(v1((ind1-1)*K+1:ind1*K));
    north=detrend(n1((ind1-1)*K+1:ind1*K));
    east=detrend(e1((ind1-1)*K+1:ind1*K));
    time=t1((ind1-1)*K+1:ind1*K);
    
    horizontalamp=zeros(fsteps,1);
    verticalamp=zeros(fsteps,1);
    horizontallist=zeros(fsteps,1);
    verticallist=zeros(fsteps,1);
    Tmax=max(time);
    thetas=zeros(fsteps,ceil(Tmax*fend));
    corr=zeros(fsteps,ceil(Tmax*fend));
    ampl=zeros(fsteps,ceil(Tmax*fend));
    dvmax=zeros(fsteps,1);

    for findex=1:1:fsteps    % loop over the frequencies
        f=fl(findex,ind1);
  
        df=dfpar*f;          % setting up the filter limits
        fmin=max(fstart,f-df/2);
        fmax=min(fnyq,f+df/2);
        flist(findex)=f;
        DT=cycles/f;
        wl=round(DT/tau);

        [na,wn]=cheb1ord([fmin+(fmax-fmin)/10,fmax-(fmax-fmin)/10]/fnyq,[fmin-(fmax-fmin)/10,fmax+(fmax-fmin)/10]/fnyq,1,5); % setting up the Chebyshev filter
        [ch1,ch2]=cheby1(na,0.5,wn);

        taper1=0:1/round(size(time,1)/100):1;taper2=ones(1,size(time,1)-size(taper1,2)*2);taper3=fliplr(taper1);
        taper=transpose([taper1,taper2,taper3]);

        % filtering the signals
        norths=filter(ch1,ch2,taper.*north);
        easts=filter(ch1,ch2,taper.*east);
        verts=filter(ch1,ch2,taper.*vert);

        derive=(sign(verts(2:K))-sign(verts(1:(K-1))))/2; % finding the negative-positive zero crossings on the vertical component

        vertsum=zeros(wl,1);
        horsum=zeros(wl,1);
        dvindex=0;

        for index=ceil(1/(4*f*tau))+1:1:length(derive)-wl % loop over all zero crossings
            if derive(index)==1
               dvindex=dvindex+1;
               vsig=verts(index:(index+wl-1));
               esig=easts(index-floor(1/(4*f*tau)):(index-floor(1/(4*f*tau))+wl-1));
               nsig=norths(index-floor(1/(4*f*tau)):(index-floor(1/(4*f*tau))+wl-1));
               integral1=sum(vsig.*esig);
               integral2=sum(vsig.*nsig);
               theta=atan(integral1/integral2);    
               if integral2<0
                  theta=theta+pi;
               end
               theta=mod(theta+pi,2*pi);  % The azimuth is well estimated in this way (assuming retrograde)
               hsig=sin(theta)*esig+cos(theta)*nsig; % The horizontal signal is projected in the azimuth direction. 
               correlation=sum(vsig.*hsig)/sqrt(sum(vsig.*vsig)*sum(hsig.*hsig)); % The correlation is always negative (between -1 and 0).
               if correlation>=-1
                  vertsum=vertsum+correlation^2*vsig;
                  horsum=horsum+correlation^2*hsig;
               end
               thetas(findex,dvindex)=theta;
               correlationlist(index)=correlation;
               thetalist(index)=theta;
               corr(findex,dvindex)=correlation;
               dvmax(findex)=dvindex;
               ampl(findex,dvindex)=sum(vsig.^2+hsig.^2);
            end
        end

        klimit=round(DT/tau);
        verticalamp(findex)=sqrt(sum(vertsum(1:klimit).^2));
        horizontalamp(findex)=sqrt(sum(horsum(1:klimit).^2));

    end
    ellist=horizontalamp./verticalamp;

    fl(:,ind1)=flist;
    el(:,ind1)=ellist;

end
V=fl;
W=el;
