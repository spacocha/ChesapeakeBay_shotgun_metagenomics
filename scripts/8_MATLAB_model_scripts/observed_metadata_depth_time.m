% Need to have run process_TPM_chesapeake_SPP

%mk_EOFs
metadata=readtable("MATLAB_metadata.txt", "ReadRowNames",true);

colormap("nebula")
clf
subplot(331)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.TEMP, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(a) Temperature (^\circ C)")
colorbar('Location','eastoutside')

subplot(332)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.SALT, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(b) Salinity (PSU)")
colorbar('Location','eastoutside')

subplot(333)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.CHLA, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(c) Chl a (mg/m^3)")
colorbar('Location','eastoutside')

subplot(334)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.DO, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(d) Oxygen (\muM)")
colorbar('Location','eastoutside')

subplot(335)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.PH, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(e) pH")
colorbar('Location','eastoutside')

subplot(336)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.PO4F, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(f) Phosphate (\mu M)")
colorbar('Location','eastoutside')

subplot(337)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.NH4F, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(g) Ammonium (\mu M)")
colorbar('Location','eastoutside')

subplot(338)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.NO3F, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(h) Nitrate (\mu M)")
colorbar('Location','eastoutside')

subplot(339)
scatter(metadata.Yearday, -1*metadata.Depth, 40, metadata.NO23F, 'filled')
xlabel('Yearday')
ylabel('Depth (m)')
title("(i) Nitrate+Nitrite (\mu M)")
colorbar('Location','eastoutside')

%Save figure from the figure GUI for best quality