%% create dataset

clear 

lon = 0:0.5:360;
lat = -81:0.5:81;
lon_n = length(lon);
lat_n = length(lat);

[xx,yy] = meshgrid(lon,lat);
yy = flipud(yy);
x = reshape(xx,lon_n*lat_n,1);
y = reshape(yy,lon_n*lat_n,1);

month = nan*zeros(lon_n*lat_n,12);

w = 2*pi/12;
w2 = 2*pi/6;

for i = 1:12
    month(:,i) = 0.2*mod(abs(y),81)*sin(w*i) + 0.1*mod(abs(x),360)*cos(w*i*2)+randn(lon_n*lat_n,1);
end

figure;
pcolor(xx,yy,reshape(month(:,1),lat_n,lon_n));

for i = 1:12
    fid = fopen(['vari.',int2str(i),'.xyz'],'w');
    fprintf(fid,'%.3f,%.3f,%6f\n',[x,y,month(:,i)]');
    fclose(fid);
end

