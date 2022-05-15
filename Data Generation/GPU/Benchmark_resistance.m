clear all;

exact = [6*pi;0;0];

points = logspace(log10(0.04),0,10);
EPS = [1e-5 3e-04 1e-2 1e-1 0.5];

storage = zeros(8,numel(points)*numel(EPS));

for i = 1:numel(points)
    disp(i+"/"+numel(points))
    gpuDevice(1);
    coarse = prolateSpheroid(1,1,points(i));
    fine = prolateSpheroid(1,1,points(i)/4);
    fine2 = [fine;coarse];
%     coarse = reshape(coarse.',[],1);
%     fine = reshape(fine.',[],1);
%     fine2 = reshape(fine2.',[],1);
    coarse = gpuArray(reshape(coarse.',[],1));
    fine = gpuArray(reshape(fine.',[],1));
    fine2 = gpuArray(reshape(fine2.',[],1));

    
    for j =1:numel(EPS)
        tic
        R = Rigid_Resistance2(coarse,coarse,[1 0 0],[0 0 0],EPS(j),1);
        errNys1 = norm(exact-R)/norm(exact);
        time1 = toc;

        tic
        R = Rigid_Resistance(coarse,fine,[1 0 0],[0 0 0],EPS(j),1);
        errNEAR = norm(exact-R)/norm(exact);
        time3 = toc;

        tic
        R = Rigid_Resistance(coarse,fine2,[1 0 0],[0 0 0],EPS(j),1);
        errNEAR2 = norm(exact-R)/norm(exact);
        time4 = toc;

        disp([points(i),EPS(j),errNys1,errNEAR,errNEAR2, time1, time3, time4])
        storage(:,(i-1)*numel(EPS) + j) = [points(i),EPS(j),errNys1,errNEAR,errNEAR2, time1, time3, time4]';
        writematrix(storage,'NearestComparion.csv')
    end
end