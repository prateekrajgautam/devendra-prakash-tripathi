clc
clear all
close all
%% user defined variables starts
% format shorte
plotcolor=['-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh';'-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh'];
%% error correcting modulator
nbhc = 12; kbhc = 4; % Codeword length and message length
hEnc = comm.BCHEncoder('CodewordLength', nbhc, 'MessageLength', kbhc);
hDec = comm.BCHDecoder('CodewordLength', nbhc, 'MessageLength', kbhc);
% msg = randi([0 1],512,1); % Four random binary messages
% 
% % Simplest syntax for encoding
% c1 = step(hEnc,msg); % BCH encoding
% d1 = step(hDec,c1); % BCH decoding
% ,
% % Check that the decoding worked correctly.
% chk = isequal(d1,msg)
%%
int=1;
papr_i=[];
papr_o=[];
papr_si=[];
papr_sl=[];
for eec=[1,2]%error correction code 1-without error correction code 2-with error correction
for q=[2]%1 for IDMA & 2 for SC-FDMA-IDMA or OFDM-IDMA
    for n1=[8]
        for it1=[3]
            for type=[1,2,3]%1 for OFDM-IDMA, 2 for SC-FDMA-IDMA(ifdma), 3 for SC-FDMA-IDMA(lfdma)
                if (int==2)&&(q==1)
                    continue
                end
                if q==1
                    typename='idma'
                end
                if q>1
                    if type==1
                        typename='OFDM-IDMA'
                    end
                    if type==2
                        typename='SC-FDMA-IDMA(ifdma)'
                    end
                    if type==3
                        typename='SC-FDMA-IDMA(lfdma)'
                    end
                end
                block=5;
                n=n1                %no.of users
                m=512;              %data length
                sl=16;              %spread length
                chiplen=m*sl;       %chip length
                itnum=it1;            %no of iteration
                ebnostart=0;        %step iteration
                ebnostep=1;
                ebnonum=10;
                useldpc=[];
                useenergyprofile=0;%set useenergyprofile to one(1) to use useenergyprofile  else set it to 2
                usemodem=1; %select modulation type '1 for psk 2 for dpsk'
                %% user defined variables ends
                hk=ones(1,n);%assume we know h_k required for idma decoder
                %% idma 
                spreadingunipolar=spreadingsequece(sl,[1,0]);%spreading sequence producing {1,0,1,0--------}
                spreading=2*spreadingunipolar-1;
                scrambrule=idmascramble(n,chiplen);  
                %% the simulation process begins
                for z=1:ebnonum
                    ebno=ebnostart+z*ebnostep;
                    snr(z)=(10.^(ebno/10))/sl;
                    sigma=sqrt(0.5/snr(z));
                    error=0;    
                    for bloc=1:block
                    %%                     transmitter
                        %% transmitter section begins
                        if eec==1
                            data=randi([0,1],n,m);
                            datatx=data;
                        end
                        if eec==2
%                             typename=cat(2,typename,'with EEC');
                            datatx=randi([0,1],n,m);
                            for e1=1:n
                                data(e1,:)=step(hEnc,datatx(e1,:)')';
                            end
                            chiplen=length(data)*sl;
                            scrambrule=idmascramble(n,chiplen); 
                        end
                        chip=spreador(spreadingunipolar,data);   
                        transmit1=interleavor(chip,scrambrule);       %transmitting data interleavor   
                        %% bpsk encoder to transmit & apply energy profile
                        t2=2*((transmit1)>0)-1;
                        %% idma multipleser (sum)
                        [r,M]=size(t2);
                        N=q*M;
                        if q==1%idma
                            tx2=sum(t2);
                        end
                        if q>1 %scfdma or ofdma
                            for i=1:r
                                X(i,:)=fft(t2(i,:));
                            end
                            if type==2
                                for i=1:r
                                    Y(i,:)=zeros(1,N);
                                    Y(i,[1:q:N])=X(i,:);%Ifdma
                                    
                                end
                                for i=1:r
                                    y1(i,:)=ifft(Y(i,:));
                                end
                                tx2=sum(y1);                          
                            end
                            if type==3
                                for i=1:r
                                    Y(i,:)=zeros(1,N);
                                    Y(i,1:length(X(i,:)))=X(i,:);%Lfdma
                                    
                                end
                                for i=1:r
                                    y1(i,:)=ifft(Y(i,:));
                                end
                                tx2=sum(y1);                          
                            end
                            if type==1
                                tx2=sum(X);
                            end
                        end
                        %% noise
                        R=channelnoise(tx2,ebno);%awgn channel
%%                         papr calc
                        if q>1%%scfdma or ofdma
                            if type==2%scfdma ifdma
    %%                          papr for SC_fdma_idma
                                Xn=abs(R);
                                Xn1=Xn.*Xn;
                                Xn2=abs(Xn1);
                                Xnmax=max(Xn2);
                                Xnavg=sum(Xn2)/numel(Xn2);
                                paprratio=Xnmax/Xnavg;
                                paprdBm=10*log10(paprratio);
                                papr_si=[papr_si,paprdBm];
                            end
                            if type==3%scfdma lfdma
    %%                          papr for SC_fdma_idma
                                Xn=abs(R);
                                Xn1=Xn.*Xn;
                                Xn2=abs(Xn1);
                                Xnmax=max(Xn2);
                                Xnavg=sum(Xn2)/numel(Xn2);
                                paprratio=Xnmax/Xnavg;
                                paprdBm=10*log10(paprratio);
                                papr_sl=[papr_sl,paprdBm];
                            end
                            
                            if type==1%ofdma
                                Xn=abs(R);
                                Xn1=Xn.*Xn;
                                Xn2=abs(Xn1);
                                Xnmax=max(Xn2);
                                Xnavg=sum(Xn2)/numel(Xn2);
                                paprratio=Xnmax/Xnavg;
                                paprdBm=10*log10(paprratio);
                                papr_o=[papr_o,paprdBm];
                            end
                        end
                        if q<2%idma
                            Xn=abs(R);
                            Xn1=Xn.*Xn;
                            Xn2=abs(Xn1);
                            Xnmax=max(Xn2);
                            Xnavg=sum(Xn2)/numel(Xn2);
                            paprratio=Xnmax/Xnavg;
                            paprdBm=10*log10(paprratio);
                            papr_i=[papr_i,paprdBm];
                        end

%                         R=tx2;%ideal channel
                        if q>1%scfdma or ofdma
                            if type==2
                                Y=fft(R);
                                X=Y(:,[1:q:N]);
                                r1=ifft(X);
                            end
                            if type==3
                                Y=fft(R);
                                X=Y(:,[1:N/q]);
                                r1=ifft(X);
                            end
                            if type==1
                                r1=ifft(R);
                            end
                        end
                        if q==1%idma
                            r1=R;
                        end
                        appllr=idmadecoderbpsk(sigma,hk,n,length(data),sl,itnum,chiplen,r1,spreading,scrambrule,useldpc);
%%                      reciever
%%                      idma decoder
                        e=0;
                        appllrf=appllr>0;% decision whether 1 or 0 is recieved
                        if eec==1
                            decodeddata=appllrf;
                        end
                        if eec==2
                            for e2=1:n
                                decodeddata(e2,:)=step(hDec,appllrf(e2,:)')';
                            end
                        end
                        [e,bertemp]=errortx(datatx,decodeddata);% check error
                        error=error+e;       
                        ber=(error/(n*m*bloc));
                        [' z=' num2str(z) ' blockno=' num2str(bloc) ' ebno=' num2str(ebno) ' ber=' num2str(ber)]
                    end
                    ebn(1,z)=ebno;
                    ber1(1,z)=ber;
                end
                result(int).ebn=ebn;
                result(int).ber=ber1;
                if q==1
                    if eec==1
                        lname{int}=[typename ',n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc)]
                    end
                    if eec==2
                        lname{int}=[typename ',n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc) 'with BHC']
                    end
                end           
                if q>1
                    if eec==1
                        lname{int}=[typename ',q=' num2str(q) ',n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc)];
                    end
                    if eec==2
                        lname{int}=[typename ',q=' num2str(q) ',n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc) 'with BHC'];
                    end
                end
                int=int+1;
                save('saveresult')
                save('temp1','q','n1','it1','int','result','plotcolor','lname','type','bloc','papr_o','papr_i','papr_si','papr_sl','hEnc','hDec','eec')
                clear all
                load temp1
            end     
        end
    end
end
end
load('saveresult')
%%  plotting result 
clc
close all
for i=1:int-1
    semilogy(result(i).ebn,(result(i).ber),plotcolor(i,:))     
    hold on
    disp(lname{i})
    disp(result(i))
end
xlabel('Eb/No')
ylabel('Bit Error Rate')
grid on
title('COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA with EEC(BHC)')
legend(lname)
saveas(gcf,'graph COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA with EEC(BHC)','jpg')
saveas(gcf,'graph COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA with EEC(BHC)','pdf')
saveas(gcf,'graph COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA with EEC(BHC)','fig')


%% plot papr
figure ()
[N1,Zi] = hist(papr_i, 10);
[N1,Zo] = hist(papr_o, 10);
[N1,Zsi] = hist(papr_si, 10);
[N1,Zsl] = hist(papr_sl, 10);
semilogy(Zi,smooth(1-cumsum(N1)/max(cumsum(N1))),'-*r')
hold on
semilogy(Zo,smooth(1-cumsum(N1)/max(cumsum(N1))),'-+g')
hold on
semilogy(Zsi,smooth(1-cumsum(N1)/max(cumsum(N1))),'-^c')
hold on
semilogy(Zsl,smooth(1-cumsum(N1)/max(cumsum(N1))),'--b')

title ('PAPR of IDMA,SC-FDMA-IDMA and OFDM-IDMA for BPSK')
xlabel ('PAPR[dB]')
ylabel ('{PAPR(PAPR>PAPR0)}')
legend(lname)
hold off
grid off;

saveas(gcf,'graph papr COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA','jpg')
saveas(gcf,'graph papr COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA','pdf')
saveas(gcf,'graph papr COMPARISION OF IDMA, SC-FDMA-IDMA & OFDM-IDMA','fig')