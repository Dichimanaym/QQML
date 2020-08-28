clear;close all
load L6to30va0to20av100f
for vi=1:41
for ei=1:81,
% figure(1);
% hold off;
% col=jet(4);
% for li=1:4;
%     semilogy(En(ei),abs(mean(Trans(ei,:,vi,li),2)),En(ei),exp(mean(log(abs(Trans(ei,:,vi,li))),2)),'color',col(li,:));
%     hold on;
% end;

%figure(2);
%hold off;
    for li=1:3
        avn(li)=mean(exp(mean(log(abs(Trans(ei,:,vi,li))),2)),1);
        avn2(li)=mean((mean((abs(Trans(ei,:,vi,li))),2)),1);
        stdn(li)=std(exp(mean(log(abs(Trans(ei,:,vi,li))),2)),1);
        stdn2(li)=std((mean((abs(Trans(ei,:,vi,li))),2)),1);
    end;
% errorbar(Ln,abs(avn),stdn);hold on
% errorbar(Ln,avn2,stdn2)
p=polyfit(Ln,abs(avn),1);
p2=polyfit(Ln,avn2,1);
gl(ei,vi)=p(1);
gl2(ei,vi)=p2(1);
end;
end;
figure(1)
%subplot(2,1,1)
surf(van,En,gl);shading interp;view(0,90);colormap(jet);colorbar;caxis([-.0001 .0001])
ylabel('E');xlabel('Disorder')
%subplot(2,1,2)
%surf(van,En,gl2);shading interp;view(0,90);colormap(jet);colorbar;caxis([-.005 .005])
%ylabel('En');xlabel('Disorder')

% load L6to30va15to20av100En0b
% Trans=gather(Trans);
% for vi=1:6
% for ei=1:41,
%     for li=1:5
%         avn(li)=mean(exp(mean(log(abs(Trans(ei,:,vi,li))),2)),1);
%         avn2(li)=mean((mean((abs(Trans(ei,:,vi,li))),2)),1);
%         stdn(li)=std(exp(mean(log(abs(Trans(ei,:,vi,li))),2)),1);
%         stdn2(li)=std((mean((abs(Trans(ei,:,vi,li))),2)),1);
%     end;
% % errorbar(Ln,abs(avn),stdn);hold on
% % errorbar(Ln,avn2,stdn2)
% p=polyfit(Ln,abs(avn),1);
% p2=polyfit(Ln,avn2,1);
% gll(ei,vi)=p(1);
% gll2(ei,vi)=p2(1);
% end;
% end;
% figure(2)
% subplot(2,1,1)
% surf(van,En,gll);shading interp;view(0,90);colormap(jet);colorbar;caxis([-.0005 .0005])
% ylabel('En');xlabel('Disorder')
% subplot(2,1,2)
% surf(van,En,smooth2b(gll2,0,0));shading interp;view(0,90);colormap(jet);colorbar;caxis([-.005 .005])
% ylabel('En');xlabel('Disorder')
% figure(3)
% plot(van,mean(gll),van,mean(gll2))
