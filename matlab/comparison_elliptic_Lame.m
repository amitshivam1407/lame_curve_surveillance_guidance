%-------------------------------------------------------------------------%
%------------------      6th March, 2024        --------------------------%
%------------------   Superellipse vector field     ----------------------%
%-------------------------------------------------------------------------%

clear;close all;clc;

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global Vg l m a_ell b_ell g kaidot_max a b n xc yc k_kai k_kaidot clockwise k_se te  kaidot_ode kaidot kaiddot_ell kaidot_ell
Vg = 15;
% l = 100;
l = 150;
m = 75;
% n = 3.6545;
% n = 3.7132;
n = 4.0631;
% a = 57.8734;
% a = 57.80;
a = 82.4549 ;
% b = 47.7350;
b = (m*a)/((2^n)*(a^n) - l^n )^(1/n) ;
xc = 0;
yc = 0;
k_kai = 100;
k_kaidot = 10;
k_se = 0.5;
clockwise = 1;

phi = deg2rad(45);
g = 9.81 ;
kaidot_max = g*tan(phi)/Vg;
% kappa_max_des = round(g*tan(phi)/Vg.^2,3) ;
kappa_max_des = 0.0436 ;
psi  = 0:2*pi/500:2*pi;

syms s
eqn = 4*(s^2) - kappa_max_des*(m^2)*s - l^2 == 0;
ss = solve(eqn,s);
a_arr = double(ss) ;
for i = 1:length(a_arr)
    if a_arr(i) > 0
        a_ell = a_arr(i) ;
        b_ell =  (m*a_ell)/((2^2)*(a_ell^2) - l^2 )^(1/2) ;
    end
end
x_ell = a_ell*cos(psi) ; y_ell = b_ell*sin(psi) ;

arclength_ellipse = get_length_ell(psi,a_ell,b_ell,0,pi/2);

xr = a*(abs(cos(psi))).^(2/n); yr = b*(abs(sin(psi))).^(2/n);
perimeter = 0;
for j = 1:length(psi)-1
    perimeter= perimeter + sqrt((xr(j) -xr(j+1))^2 + (yr(j) -yr(j+1))^2);
end
% ae_arr(i) = x(i,1);
arclength_se = perimeter;

arclength_percentage = ((arclength_ellipse - arclength_se)./arclength_ellipse)*100;


psi  = 0:2*pi/500:2*pi;
range = -160:5:160;
[X,Y] = meshgrid(range);

tspan = [0 50];  % run time of ode
%----------- initial value ----------
x0   = -160;
y0   = -140;
[kai0,kaidot0] = get_initialvalue(x0,y0);
[kaid0_ell,kaidot0_ell]= get_intialellipse(x0,y0);
x_init = [x0;y0;kai0;kaidot0;x0;y0;kaid0_ell;kaidot0_ell];
% x_init = [x0;y0;kai0];
%%---------- Ode function ----------------------------------------------%

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,x]   = ode45(@(t,x) fun_ellipse_se(t,x), tspan, x_init,options);

x_ini     = x(1,1);
y_ini     = x(1,2);
x_end     = x(end,1);
y_end     = x(end,2);

Tracking_error = fun_error(x(:,1),x(:,2));
xdot = Vg*cos(x(:,3));
ydot = Vg*sin(x(:,3));
[kaid,kaid_dot] = get_kai(x(:,1),x(:,2),xdot,ydot);

% plot parameter
ax_fnt = 18;
lbl_fnt = 21;
ax_wdth = 3;
lgd_fnt = 17;

figure(1)
[x_rect, y_rect] = get_rectangle(l,m);
plot(x_rect,y_rect,'-','linewidth',3,'MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','magenta'); hold on;
% rectangular_boundary = get_rectangle(l,m);
superellipse_path = circumscribed_se(psi);hold on;
% [xdot_out,ydot_out,X_out,Y_out] = vf_superellipse(X,Y);
[X_out,Y_out,X_in,Y_in,xdot_out,ydot_out,xdot_in,ydot_in] = get_superellipsevf(X,Y);
quiver(X_out,Y_out,xdot_out,ydot_out,'Color',[0.75 0.75 0.75],'LineWidth',1);hold on;
% quiver(X_in,Y_in,xdot_in,ydot_in,'Color',[0.75 0.75 0.75],'LineWidth',1);hold on;
% h13 = plot(xc,yc,'--rs','LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',10);hold on;
% h13.Annotation.LegendInformation.IconDisplayStyle = 'off';
h1 = plot(x(:,1),x(:,2),'r','linewidth',3);hold on ;
h15 = plot(x(:,5),x(:,6),'b','linewidth',3);hold on ;
h11 = plot(x_ini,y_ini,'-o','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','green'); hold on;
h11.Annotation.LegendInformation.IconDisplayStyle = 'off';
h12=plot(x_end,y_end,'-s','linewidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','cyan'); hold on;
h12.Annotation.LegendInformation.IconDisplayStyle = 'off';

ax1 = gca;
ax1.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',3)
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
legend(ax1,'Rectangle','Circumscribed path','Vector field','UAV trajectory','Fontsize',lgd_fnt,'NumColumns',1)
axis(ax1, 'equal');

figure(2)
h2 = plot(t,Tracking_error,'linewidth',3);hold on ;grid on;
ax2 = gca;
ax2.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',3)
xlabel(ax2,'  $t$, s','Fontsize',lbl_fnt);
ylabel(ax2,'  $(\beta - 1) $, m','Fontsize',lbl_fnt);


figure(3)
h31 = plot(t,wrapToPi(kaid)*(180/pi),'r','linewidth',3);hold on ;grid on;
h32 = plot(t,wrapToPi(x(:,3))*(180/pi),'--b','linewidth',3);hold on ;grid on;
ax3 = gca;
ax3.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',3)
xlabel(ax3,' Time $t$, s','Fontsize',lbl_fnt);
ylabel(ax3,' Course angle $\chi $, deg.','Fontsize',lbl_fnt);
legend(ax3,'Commanded','Achieved','Fontsize',lgd_fnt)

figure(4)
h41 = yline(kaidot_max,'--k','linewidth',3);hold on ;grid on;
h44 = yline(-kaidot_max,'--k','linewidth',3);hold on ;grid on;
h44.Annotation.LegendInformation.IconDisplayStyle = 'off';
h42 = plot(te,kaidot_ode,'r','linewidth',3);hold on ;grid on;
h43 = plot(te,kaidot,'--b','linewidth',3);hold on ;grid on;
% h42 = plot(t,kaid_dot,'r','linewidth',3);hold on ;grid on;
% h4.Color = blue;
% h43 = plot(t,x(:,4),'b','linewidth',3);hold on ;grid on;
ax4 = gca;
ax4.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax4.XColor = 'black';         % Box horizontal lines' color
ax4.YColor = 'black';         % Box vertical lines' color
set(ax4,'linewidth',3)
xlabel(ax4,' $t$, s','Fontsize',lbl_fnt);
yticks(ax4,[-1.25 -0.6533 0 0.6533 1.25])
yticklabels(ax4,{'-1.25','-0.6533','0','0.6533','1.25'})
ylabel(ax4,' $\dot{\psi} $, rad/s','Fontsize',lbl_fnt);
legend(ax4,' $\dot{\psi}_{\mathrm{max}}$','Commanded','Achieved','Fontsize',lgd_fnt)
% legend(ax4,'Commanded','Achieved','Fontsize',lgd_fnt)
axis(ax4,[0 t(end) -kaidot_max-0.5 kaidot_max+0.5])

figure(5)
h51 = plot(x_rect,y_rect,'-d','linewidth',3,'DisplayName','Rectangle'); hold on;
% h51.Annotation.LegendInformation.IconDisplayStyle = 'on';
% h52 = plot(xc,yc,'-o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor','green','DisplayName','Center '); hold on;
% h52.Annotation.LegendInformation.IconDisplayStyle = 'off';
superellipse_path = circumscribed_se(psi);hold on;
h53 = plot(x_ell,y_ell,'-k','linewidth',3,'DisplayName','Circumscribing ellipse'); hold on;
% h53.Annotation.LegendInformation.IconDisplayStyle = 'on';
ax5 = gca;
ax5.FontSize = ax_fnt;
box on                      % Switch on the box around the axis
ax5.XColor = 'black';         % Box horizontal lines' color
ax5.YColor = 'black';         % Box vertical lines' color
set(ax5,'linewidth',ax_wdth) ;
xlabel(ax5,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax5,'$ y, $ m','Fontsize',lbl_fnt);
legend(ax5,'Rectangle',"Lam\'e curve",'Ellipse','Fontsize',lgd_fnt)
axis(ax5,'equal')

figure(6)
h61 = yline(kaidot_max,'--k','linewidth',3);hold on ;grid on;
h62 = yline(-kaidot_max,'--k','linewidth',3);hold on ;grid on;
h62.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h42 = plot(te,kaidot_ode,'r','linewidth',3);hold on ;grid on;
h63 = plot(te,kaidot,'r','linewidth',3);hold on ;grid on;
h64 = plot(te,kaidot_ell,'b','linewidth',3);hold on ;grid on;
% h61 = plot(t,kaid_dot./Vg,'r','linewidth',3);hold on ;grid on;
% h62 = plot(te,-kaiddotell./Vg,'b','linewidth',3);hold on ;grid on;
% h62 = plot(t,x(:,8)./Vg,'b','linewidth',3);hold on ;grid on;
ax6 = gca;
ax6.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax6.XColor = 'black';         % Box horizontal lines' color
ax6.YColor = 'black';         % Box vertical lines' color
set(ax6,'linewidth',3)
xlabel(ax6,'  $t$, s','Fontsize',lbl_fnt);
yticks(ax6,[-1.25 -0.6533 0 0.6533 1.25])
yticklabels(ax6,{'-1.25','-0.6533','0','0.6533','1.25'})
ylabel(ax6,' $\dot{\psi} $, rad/s','Fontsize',lbl_fnt);
legend(ax6,' $\dot{\psi}_{\mathrm{max}}$',"Lam\'e curve",'Ellipse','Fontsize',lgd_fnt)
% axis(ax6,[0 t(end) -kaidot_max-0.5 kaidot_max+0.5])
axis(ax6,[0 t(end) -1.3 1.3])

%% Create GIF for Figure(6): Turn Rate Comparison
gif_name_f6 = 'turnrate_comparison.gif';
% Delete existing GIF file to avoid format compatibility issues
if exist(gif_name_f6, 'file')
    delete(gif_name_f6);
end

figure(60); clf;
set(gcf,'Color','w','Position',[100 100 720 400]);
ax_gif6 = gca;
hold(ax_gif6,'on');
grid(ax_gif6,'on');
box(ax_gif6,'on');

% Static background: limit lines
yline(ax_gif6, kaidot_max, '--k', 'LineWidth', 3);
yline(ax_gif6, -kaidot_max, '--k', 'LineWidth', 3);

% Animated curves (start empty)
h_se_dot = plot(ax_gif6, nan, nan, 'r', 'LineWidth', 3);
h_ell_dot = plot(ax_gif6, nan, nan, 'b', 'LineWidth', 3);

xlabel(ax_gif6,'  $t$, s','Fontsize',lbl_fnt);
ylabel(ax_gif6,' $\dot{\psi} $, rad/s','Fontsize',lbl_fnt);
title(ax_gif6,'Turn Rate Comparison','Fontsize',lbl_fnt);

ax_gif6.FontSize = ax_fnt;
ax_gif6.XColor = 'black';
ax_gif6.YColor = 'black';
set(ax_gif6,'linewidth',ax_wdth)
axis(ax_gif6,[0 t(end) -1.3 1.3]);

legend(ax_gif6,' $\dot{\psi}_{\mathrm{max}}$','',"Lam\'e curve",'Ellipse','Fontsize',lgd_fnt,'Location','northeast');

% choose fewer frames for smaller GIF
skip6 = 10;
delayTime6 = 0.08;

for k = 1:skip6:length(te)
    set(h_se_dot, 'XData', te(1:k), 'YData', kaidot(1:k));
    set(h_ell_dot, 'XData', te(1:k), 'YData', kaidot_ell(1:k));
    title(ax_gif6, sprintf('Time: %.2f s', te(k)), 'Fontsize', lbl_fnt);
    drawnow;

    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,64);

    if k == 1
        imwrite(A,map,gif_name_f6,'gif','LoopCount',Inf,'DelayTime',delayTime6);
    else
        imwrite(A,map,gif_name_f6,'gif','WriteMode','append','DelayTime',delayTime6);
    end
end

fprintf('Turn rate GIF saved to: %s\n', gif_name_f6);

%% Figure 7: Animation with Fixed-Wing UAV Markers
figure(7)
fig = gcf;
set(fig, 'Position', [100 100 900 650]);
ax7 = gca;
[x_rect, y_rect] = get_rectangle(l,m);
plot(x_rect,y_rect,'-','linewidth',3); hold on; grid on;
psi_anim = 0:2*pi/500:2*pi;
xd_anim = xc + a*(abs(cos(psi_anim)).^(2/n)).*sign(cos(psi_anim));
yd_anim = yc + b*(abs(sin(psi_anim)).^(2/n)).*sign(sin(psi_anim));
plot(xd_anim,yd_anim,'g','LineWidth',3); hold on;
plot(x_ell,y_ell,'-k','linewidth',3); hold on;
quiver(X_out,Y_out,xdot_out,ydot_out,'Color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(xc,yc,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','magenta'); hold on;
plot(x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10); hold on;
ax7.FontSize = ax_fnt;
box on
ax7.XColor = 'black';
ax7.YColor = 'black';
set(ax7,'linewidth',ax_wdth)
xlabel(ax7,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax7,'$ y, $ m','Fontsize',lbl_fnt);
axis(ax7,'equal');

% Initialize trajectory trails and UAV markers
h_trail1 = plot(nan, nan, 'r', 'linewidth', 3);
h_trail2 = plot(nan, nan, 'b', 'linewidth', 3);

% Create fixed-wing UAV shapes
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav1 = fill(nan, nan, 'g', 'EdgeColor', 'k', 'LineWidth', 2);
h_uav2 = fill(nan, nan, 'c', 'EdgeColor', 'k', 'LineWidth', 2);

legend('Boundary',"Lam\'e curve circumscription",'Elliptic circumscription','Vector field','','',"Lam\'e curve trajectory",'Elliptic trajectory','','','Fontsize',lgd_fnt,'NumColumns',2);

% Setup video writer
video_filename = sprintf('ellipse_superellipse_animation_%s.mp4', datestr(now,''));
vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 10;
vidObj.Quality = 95;
open(vidObj);

% Animation loop
dt_anim = 0.1;
t_anim = 0:dt_anim:t(end);
for i = 1:length(t_anim)
    % Interpolate positions and headings
    x1_curr = interp1(t, x(:,1), t_anim(i));
    y1_curr = interp1(t, x(:,2), t_anim(i));
    chi1_curr = interp1(t, x(:,3), t_anim(i));
    x2_curr = interp1(t, x(:,5), t_anim(i));
    y2_curr = interp1(t, x(:,6), t_anim(i));
    chi2_curr = interp1(t, x(:,7), t_anim(i));
    
    % Update trajectory trails
    idx = find(t <= t_anim(i));
    set(h_trail1, 'XData', x(idx,1), 'YData', x(idx,2));
    set(h_trail2, 'XData', x(idx,5), 'YData', x(idx,6));
    
    % Rotate and translate UAV shapes
    R1 = [cos(chi1_curr), -sin(chi1_curr); sin(chi1_curr), cos(chi1_curr)];
    uav1_rotated = R1 * [uav_shape_x; uav_shape_y];
    uav1_x = uav1_rotated(1,:) + x1_curr;
    uav1_y = uav1_rotated(2,:) + y1_curr;
    set(h_uav1, 'XData', uav1_x, 'YData', uav1_y);
    
    R2 = [cos(chi2_curr), -sin(chi2_curr); sin(chi2_curr), cos(chi2_curr)];
    uav2_rotated = R2 * [uav_shape_x; uav_shape_y];
    uav2_x = uav2_rotated(1,:) + x2_curr;
    uav2_y = uav2_rotated(2,:) + y2_curr;
    set(h_uav2, 'XData', uav2_x, 'YData', uav2_y);
    
    title(ax7, sprintf('Time: %.2f s', t_anim(i)), 'Fontsize', lbl_fnt);
    drawnow;
    
    frame = getframe(gcf);
    writeVideo(vidObj, frame);
end

close(vidObj);
fprintf('Animation saved to: %s\n', video_filename);
disp('Animation complete!');

%% Create GIF animation: UAVs moving on vector field
gif_name = 'ellipse_superellipse_path_following.gif';
% Delete existing GIF file to avoid format compatibility issues
if exist(gif_name, 'file')
    delete(gif_name);
end

figure(10); clf;
set(gcf,'Color','w','Position',[100 100 720 480]);
ax_gif = gca;
hold(ax_gif,'on');
grid(ax_gif,'on');
box(ax_gif,'on');

% Plot fixed background once
[x_rect, y_rect] = get_rectangle(l,m);
plot(ax_gif,x_rect,y_rect,'-','linewidth',3); hold(ax_gif,'on');
psi_gif = 0:2*pi/500:2*pi;
xd_gif = xc + a*(abs(cos(psi_gif)).^(2/n)).*sign(cos(psi_gif));
yd_gif = yc + b*(abs(sin(psi_gif)).^(2/n)).*sign(sin(psi_gif));
plot(ax_gif,xd_gif,yd_gif,'g','linewidth',3); hold(ax_gif,'on');
plot(ax_gif,x_ell,y_ell,'-k','linewidth',3); hold(ax_gif,'on');
quiver(ax_gif,X_out,Y_out,xdot_out,ydot_out,'Color',[0.75 0.75 0.75],'LineWidth',1); hold(ax_gif,'on');
plot(ax_gif,xc,yc,'-o','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','magenta'); hold(ax_gif,'on');
plot(ax_gif,x_ini,y_ini,'-o','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10); hold(ax_gif,'on');

xlabel(ax_gif,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax_gif,' $ y, $ m','Fontsize',lbl_fnt);
% title(ax_gif,'Superellipse and Ellipse Path Following','Fontsize',lbl_fnt);

ax_gif.FontSize = ax_fnt;
ax_gif.XColor = 'black';
ax_gif.YColor = 'black';
set(ax_gif,'linewidth',ax_wdth)
axis(ax_gif,'equal');
% xlim(ax_gif,[-160 160]);
% ylim(ax_gif,[-160 160]);

% animated objects
traj_line1 = plot(ax_gif,nan,nan,'r','linewidth',3);
traj_line2 = plot(ax_gif,nan,nan,'b','linewidth',3);

% Create fixed-wing UAV shapes
uav_size = 5;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];
h_uav1_gif = fill(ax_gif,nan, nan, 'g', 'EdgeColor', 'k', 'LineWidth', 2);
h_uav2_gif = fill(ax_gif,nan, nan, 'c', 'EdgeColor', 'k', 'LineWidth', 2);

% legend(ax_gif,'Boundary','Lam''e curve path','Ellipse','Vector field','','','UAV-SE','UAV-Ell','','', ...
%     'Fontsize',lgd_fnt,'Location','northwest');

legend(ax_gif,'Boundary',"Lam\'e curve circumscription",'Elliptic circumscription','Vector field','','',"Lam\'e curve following",'Elliptic path following','','','Fontsize',lgd_fnt,'NumColumns',2);


% choose fewer frames for smaller GIF
skip = 10;
delayTime = 0.08;

for k = 1:skip:length(t)
    
    set(traj_line1,'XData',x(1:k,1),'YData',x(1:k,2));
    set(traj_line2,'XData',x(1:k,5),'YData',x(1:k,6));
    
    % Get current headings
    chi1_curr = x(k,3);
    chi2_curr = x(k,7);
    
    % Rotate and translate UAV shapes
    R1 = [cos(chi1_curr), -sin(chi1_curr); sin(chi1_curr), cos(chi1_curr)];
    uav1_rotated = R1 * [uav_shape_x; uav_shape_y];
    uav1_x = uav1_rotated(1,:) + x(k,1);
    uav1_y = uav1_rotated(2,:) + x(k,2);
    set(h_uav1_gif, 'XData', uav1_x, 'YData', uav1_y);
    
    R2 = [cos(chi2_curr), -sin(chi2_curr); sin(chi2_curr), cos(chi2_curr)];
    uav2_rotated = R2 * [uav_shape_x; uav_shape_y];
    uav2_x = uav2_rotated(1,:) + x(k,5);
    uav2_y = uav2_rotated(2,:) + x(k,6);
    set(h_uav2_gif, 'XData', uav2_x, 'YData', uav2_y);
    
    title(ax_gif, sprintf('Time: %.2f s', t(k)), 'Fontsize', lbl_fnt);
    drawnow;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,64);
    
    if k == 1
        imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',delayTime);
    else
        imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',delayTime);
    end
end

fprintf('GIF animation saved to: %s\n', gif_name);

function out = fun_ellipse_se(t,x)
global Vg a  b n a_ell b_ell kaidot_max k_kai k_kaidot k_se  te kaidot_ode kaidot kaiddot_ell kaidot_ell

if t==0
    te = [];
    kaidot_ode = [];
    kaidot = [];
    kaiddot_ell = [];
    kaidot_ell = [];
end

out(1,1) = Vg*cos(x(3));
out(2,1) = Vg*sin(x(3));


if (x(1) > 0) && (x(2) > 0)
    beta = (x(1).^n) /(a^n) + (x(2).^n) /(b^n) ;
    fx = (n*(x(1)^(n-1)))./(a^n) ;
    fy = (n*(x(2)^(n-1)))./(b^n) ;
    fxx = (n*(n-1)*(x(1)^(n-2)))./(a^n) ;
    fyy = (n*(n-1)*(x(2)^(n-2)))./(b^n) ;
    kait = pi - atan2(abs(fx),abs(fy)) ;
    kait_dot = -((fy*fxx*out(1,1) - fx*fyy*out(2,1))./fy^2)*(cos(kait))^2 ;
    beta_dot = fx*out(1,1) + fy*out(2,1) ;
    kaio_dot = ((2*k_se)./((1 + k_se*(beta - 1)^2)*sqrt((k_se*(beta - 1))^2 + 2*k_se)))*beta_dot;
elseif (x(1) < 0) && (x(2) > 0)
    beta = ((-x(1)).^n) /(a^n) + (x(2).^n) /(b^n) ;
    fx = -(n*((-x(1))^(n-1)))./(a^n) ;
    fy = (n*(x(2)^(n-1)))./(b^n) ;
    fxx = (n*(n-1)*((-x(1))^(n-2)))./(a^n) ;
    fyy = (n*(n-1)*(x(2)^(n-2)))./(b^n) ;
    kait = pi + atan2(abs(fx),abs(fy)) ;
    kait_dot = -((fy*fxx*out(1,1) - fx*fyy*out(2,1))./fy^2)*(cos(kait))^2 ;
    beta_dot = fx*out(1,1) + fy*out(2,1) ;
    kaio_dot = ((2*k_se)./((1 + k_se*(beta - 1)^2)*sqrt((k_se*(beta - 1))^2 + 2*k_se)))*beta_dot;
elseif (x(1) < 0) && (x(2) < 0)
    beta = ((-x(1)).^n) /(a^n) + ((-x(2)).^n) /(b^n) ;
    fx = -(n*((-x(1))^(n-1)))./(a^n) ;
    fy = -(n*((-x(2))^(n-1)))./(b^n) ;
    fxx = (n*(n-1)*((-x(1))^(n-2)))./(a^n) ;
    fyy = (n*(n-1)*((-x(2))^(n-2)))./(b^n) ;
    kait = 2*pi - atan2(abs(fx),abs(fy)) ;
    kait_dot = -((fy*fxx*out(1,1) - fx*fyy*out(2,1))./fy^2)*(cos(kait))^2 ;
    beta_dot = fx*out(1,1) + fy*out(2,1) ;
    kaio_dot = ((2*k_se)./((1 + k_se*(beta - 1)^2)*sqrt((k_se*(beta - 1))^2 + 2*k_se)))*beta_dot;
elseif (x(1) > 0) && (x(2) < 0)
    beta = (x(1).^n) /(a^n) + ((-x(2)).^n) /(b^n) ;
    fx = (n*(x(1)^(n-1)))./(a^n) ;
    fy = -(n*((-x(2))^(n-1)))./(b^n) ;
    fxx = (n*(n-1)*(x(1)^(n-2)))./(a^n) ;
    fyy = (n*(n-1)*((-x(2))^(n-2)))./(b^n) ;
    kait =  atan2(abs(fx),abs(fy)) ;
    kait_dot = -((fy*fxx*out(1,1) - fx*fyy*out(2,1))./fy^2)*(cos(kait))^2 ;
    beta_dot = fx*out(1,1) + fy*out(2,1) ;
    kaio_dot = ((2*k_se)./((1 + k_se*(beta - 1)^2)*sqrt((k_se*(beta - 1))^2 + 2*k_se)))*beta_dot;
end
% 
kai_o = pi/2 -  asin(1./(1 + k_se*(beta-1).^2))  ;


if beta <=1
     kaid = wrapToPi(kait  -  kai_o)  ;
    kaid_dot = kait_dot - kaio_dot;
else
     kaid = wrapToPi(kait  +  kai_o)  ;
    kaid_dot = kait_dot + kaio_dot;
end

if (kaid_dot < -kaidot_max)
    kaid_dot = -kaidot_max;
    
end
if (kaid_dot > kaidot_max) 
    kaid_dot = kaidot_max;    
end


out(3,1) = x(4);

if x(4) < -kaidot_max 
    x(4) = -kaidot_max ;
end
if x(4) > kaidot_max
    x(4) = kaidot_max ;
end

out(4,1) =  k_kai*wrapToPi(kaid - x(3))  + k_kaidot*(kaid_dot - x(4))  ;




beta_ell = x(5)^2/a_ell^2 + x(6)^2/b_ell^2 ;
fx_e = 2*x(5)/a_ell^2 ;
fxx_e = 2./a_ell^2 ;
fy_e = 2*x(6)/b_ell^2 ;
fyy_e = 2/b_ell^2 ;
kait_e =  atan2((b_ell^2*x(5) ),(-a_ell^2*x(6)));

kaio_e = pi/2 -  asin(1./(1 + k_se*(beta_ell-1).^2))  ;
if beta_ell <=1
    kaid_ell = kait_e - kaio_e ;
    out(5,1) = Vg*cos(x(7));
    out(6,1) = Vg*sin(x(7));
    kaitdot_ell = -((fy_e*fxx_e*out(5,1) - fx_e*fyy_e*out(6,1))./fy_e^2)*(cos(kait_e))^2 ;
    betadot_ell = fx_e*out(5,1) + fy_e*out(6,1) ;
    kaiodot_ell = ((2*k_se)./((1 + k_se*(beta_ell - 1)^2)*sqrt((k_se*(beta_ell - 1))^2 + 2*k_se)))*betadot_ell;
    kaiddotell = kaitdot_ell - kaiodot_ell;
else
    kaid_ell = kait_e + kaio_e ;
    out(5,1) = Vg*cos(x(7));
    out(6,1) = Vg*sin(x(7));
    kaitdot_ell = -((fy_e*fxx_e*out(5,1) - fx_e*fyy_e*out(6,1))./fy_e^2)*(cos(kait_e))^2 ;
    betadot_ell = fx_e*out(5,1) + fy_e*out(6,1) ;
    kaiodot_ell = ((2*k_se)./((1 + k_se*(beta_ell - 1)^2)*sqrt((k_se*(beta_ell - 1))^2 + 2*k_se)))*betadot_ell;
    kaiddotell = kaitdot_ell + kaiodot_ell;
end

if (kaiddotell < -kaidot_max)
    kaiddotell = -kaidot_max;    
end
if (kaiddotell > kaidot_max) 
    kaiddotell = kaidot_max;    
end

out(7,1) = x(8);
if x(8) < -kaidot_max 
    x(8) = -kaidot_max ;
end
if x(8) > kaidot_max
    x(8) = kaidot_max ;
end
out(8,1) = k_kai*wrapToPi(kaid_ell - x(7)) + k_kaidot*(kaiddotell - x(8)) ;




te = [te t];
kaidot_ode = [kaidot_ode kaid_dot];
kaidot = [kaidot x(4)];

kaiddot_ell = [kaiddot_ell kaiddotell];
kaidot_ell = [kaidot_ell x(8)];
end


function [x_rect,y_rect] = get_rectangle(l,m)
global xc yc 
%% Rectangular boundary --------------------------------------------------%
x1 = xc - l/2 ; y1 = yc - m/2 ; x2 = xc + l/2 ; y2 = yc -m/2 ;
x3 = xc + l/2 ; y3 = yc + m/2 ; x4 = xc - l/2 ; y4 = yc +m/2 ;
x_rect = [x1 x2 x3 x4 x1];
y_rect = [y1 y2 y3 y4 y1];
% rectangular_boundary = plot(x_rect,y_rect,'m','LineWidth',3);
end

function arc_length = get_length_ell(theta,aell,bell,min,max)
% global l m a b

% a = aa;
% b = (m*a)/((2^2)*(a^2) - l^2 )^(1/2) ;
% b = k.*a;
fun = @(theta) sqrt((aell.*cos(theta)).^2 + (bell.*sin(theta)).^2) ;
arc_length = 4*integral(fun,min,max) ;

end

function superellipse_path = circumscribed_se(psi)
global a b n xc yc
xd = xc + a*(abs(cos(psi)).^(2/n)).*sign(cos(psi)) ;
yd = yc + b*(abs(sin(psi)).^(2/n)).*sign(sin(psi)) ;
superellipse_path = plot(xd,yd,'r','LineWidth',3);
end


function [X_out,Y_out,X_in,Y_in,xdot_out,ydot_out,xdot_in,ydot_in] = get_superellipsevf(X,Y)
global Vg a b n k_se clockwise
for i = 1: length(X)
    for j = 1:length(X)
        %--------------------   vector field  ----------------------------%
        if clockwise ==0
            if X(i,j)>=0  && Y(i,j)>=0
                beta_prop(i,j) = X(i,j).^n /a^n + Y(i,j).^n /b^n ;
                des_angle(i,j) = -atan2((b^n*X(i,j).^(n-1) ),(a^n*Y(i,j).^(n-1)));
            elseif  X(i,j)<0 && Y(i,j)>=0
                beta_prop(i,j) = (-X(i,j)).^n /a^n + Y(i,j).^n /b^n ;
                des_angle(i,j) =   atan2((b^n*(-X(i,j)).^(n-1)),(a^n*Y(i,j).^(n-1)));
            elseif  X(i,j)<=0 &&  Y(i,j)<0
                beta_prop(i,j) = (-X(i,j)).^n /a^n + (-Y(i,j)).^n /b^n ;
                des_angle(i,j) = pi - atan2((b^n*(-X(i,j)).^(n-1) ),(a^n*(-Y(i,j)).^(n-1)));
            elseif X(i,j)>=0  &&  Y(i,j)<0
                beta_prop(i,j) = X(i,j).^n /a^n + (-Y(i,j)).^n /b^n ;
                des_angle(i,j) =  -pi + atan2((b^n*X(i,j).^(n-1) ),(a^n*(-Y(i,j)).^(n-1)));
            end
            kai_o(i,j) = pi/2 -  asin(1./(1 + k_se*(beta_prop(i,j)-1).^2))  ;
            if beta_prop(i,j) <1
                kaid_in(i,j) = des_angle(i,j)  +  kai_o(i,j)  ;
                X_in(i,j) = X(i,j);
                Y_in(i,j) = Y(i,j);
                xdot_in(i,j) = Vg*cos(kaid_in(i,j));
                ydot_in(i,j) = Vg*sin(kaid_in(i,j));
            else
                kaid_out(i,j) = des_angle(i,j)  -  kai_o(i,j)  ;
                X_out(i,j) = X(i,j);
                Y_out(i,j) = Y(i,j);
                xdot_out(i,j) = Vg*cos(kaid_out(i,j));
                ydot_out(i,j) = Vg*sin(kaid_out(i,j));
            end
        else
            if (X(i,j)>=0 && Y(i,j)>=0)
                beta_prop(i,j) = X(i,j).^n /a^n + Y(i,j).^n /b^n ;
%                 des_angle(i,j) = (pi + atan2(-(b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1))));
                 des_angle(i,j) = (pi - atan2((b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1))));
            elseif (X(i,j)<=0 && Y(i,j)>=0)
                beta_prop(i,j) = (-X(i,j)).^n /a^n + Y(i,j).^n /b^n ;
%                 des_angle(i,j) =  -pi - atan2(-(b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1)));
%                  des_angle(i,j) =  pi + atan2((b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1)));
                  des_angle(i,j) =  pi + atan2((b^n*(-X(i,j)).^(n-1)),(a^n*(Y(i,j)).^(n-1)));
            elseif (X(i,j)<=0 && Y(i,j)<=0)
                beta_prop(i,j) = (-X(i,j)).^n /a^n + (-Y(i,j)).^n /b^n ;
%                 des_angle(i,j) =   atan2(-(b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1)));
%                  des_angle(i,j) =  2*pi -atan2((b^n*abs(X(i,j)).^(n-1)),(a^n*abs(Y(i,j)).^(n-1)));
                 des_angle(i,j) =  2*pi -atan2((b^n*(-X(i,j)).^(n-1)),(a^n*(-Y(i,j)).^(n-1)));
            elseif (X(i,j)>=0 && Y(i,j)<=0)
                beta_prop(i,j) = X(i,j).^n /a^n + (-Y(i,j)).^n /b^n ;
%                 des_angle(i,j) =   -atan2(-(b^n*abs(X(i,j)).^(n-1) ),(a^n*abs(Y(i,j)).^(n-1)));
%                 des_angle(i,j) =   atan2((b^n*abs(X(i,j)).^(n-1) ),(a^n*abs(Y(i,j)).^(n-1)));
                des_angle(i,j) =   atan2((b^n*(X(i,j)).^(n-1) ),(a^n*(-Y(i,j)).^(n-1)));
            end
            kai_o(i,j) = pi/2 -  asin(1./(1 + k_se*(beta_prop(i,j)-1).^2))  ;
            if beta_prop(i,j) <=1
                kaid_in(i,j) = des_angle(i,j)  -  kai_o(i,j)  ;
%                 kaid_in(i,j) = des_angle(i,j)    ;
                X_in(i,j) = X(i,j);
                Y_in(i,j) = Y(i,j);
                xdot_in(i,j) = Vg*cos(kaid_in(i,j));
                ydot_in(i,j) = Vg*sin(kaid_in(i,j));
            else
                 kaid_out(i,j) = des_angle(i,j)  +  kai_o(i,j)  ;
%                 kaid_out(i,j) = des_angle(i,j)    ;
                X_out(i,j) = X(i,j);
                Y_out(i,j) = Y(i,j);
                xdot_out(i,j) = Vg*cos(kaid_out(i,j));
                ydot_out(i,j) = Vg*sin(kaid_out(i,j));
            end

        end
    end
end
end

function [kaid,kaid_dot] = get_kai(x,y,xdot,ydot)
global  a b n k_se clockwise xc yc

%------------------------ initial kai  ------------------------------%
for i = 1: length(x)
%     for j = 1:length(X)
        %--------------------   vector field  ----------------------------%
        if clockwise ==0
            if (x(i)>=0 && y(i)>=0)
                beta(i) = x(i).^n /a^n + y(i).^n /b^n ;
                fx(i) = (n*(x(i)^(n-1)))./(a^n) ;
                fy(i) = (n*(y(i)^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*(x(i)^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*(y(i)^(n-2)))./(b^n) ;
                % kai_t = (pi + atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));                
                % kai_t(i) = (pi - atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
                kai_t(i) =   pi -atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)<=0 && y(i)>=0)
                beta(i) = (-x(i)).^n /a^n + y(i).^n /b^n ;
                fx(i) = -(n*((-x(i))^(n-1)))./(a^n) ;
                fy(i) = (n*(y(i)^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*((-x(i))^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*(y(i)^(n-2)))./(b^n) ;
                % kai_t =  -pi - atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))); 
                % kai_t =  pi + atan2((b^n*(-x(i)).^(n-1)),(a^n*(y(i)).^(n-1)));               
                % kai_t(i) =  pi + atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   pi + atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)<=0 && y(i)<=0)
                beta(i) = (-x(i)).^n /a^n + (-y(i)).^n /b^n ;
                fx(i) = -(n*((-x(i))^(n-1)))./(a^n) ;
                fy(i) = -(n*((-y(i))^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*((-x(i))^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*((-y(i))^(n-2)))./(b^n) ;
                % kai_t =   atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                % kai_t =  2*pi -atan2((b^n*(-x(i)).^(n-1)),(a^n*(-y(i)).^(n-1)));               
                % kai_t(i) =  2*pi -atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   2*pi -atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)>=0 && y(i)<=0)
                beta(i) = x(i).^n /a^n + (-y(i)).^n /b^n ;
                fx(i) = (n*(x(i)^(n-1)))./(a^n) ;
                fy(i) = -(n*((-y(i))^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*(x(i)^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*((-y(i))^(n-2)))./(b^n) ;
                % kai_t =   -atan2(-(b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));                
                % kai_t =   atan2((b^n*(x(i)).^(n-1) ),(a^n*(-y(i)).^(n-1)));           
                % kai_t(i) =   atan2((b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            end
            kai_o(i) = pi/2 -  asin(1./(1 + k_se*(beta(i)-1).^2))  ;


            if beta <=1
                 kaid(i) = kai_t(i)  -  kai_o(i)  ;
                 kaid_dot(i) = kait_dot(i) - kaio_dot(i) ;
            else
                 kaid(i) = kai_t(i)  +  kai_o(i)  ;
                 kaid_dot(i) = kait_dot(i) + kaio_dot(i) ;

            end

        else
              if (x(i)>=0 && y(i)>=0)
                beta(i) = x(i).^n /a^n + y(i).^n /b^n ;
                fx(i) = (n*(x(i)^(n-1)))./(a^n) ;
                fy(i) = (n*(y(i)^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*(x(i)^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*(y(i)^(n-2)))./(b^n) ;
                % kai_t = (pi + atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));                
                % kai_t(i) = (pi - atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
                kai_t(i) =   pi -atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)<=0 && y(i)>=0)
                beta(i) = (-x(i)).^n /a^n + y(i).^n /b^n ;
                fx(i) = -(n*((-x(i))^(n-1)))./(a^n) ;
                fy(i) = (n*(y(i)^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*((-x(i))^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*(y(i)^(n-2)))./(b^n) ;
                % kai_t =  -pi - atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))); 
                % kai_t =  pi + atan2((b^n*(-x(i)).^(n-1)),(a^n*(y(i)).^(n-1)));               
                % kai_t(i) =  pi + atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   pi + atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)<=0 && y(i)<=0)
                beta(i) = (-x(i)).^n /a^n + (-y(i)).^n /b^n ;
                fx(i) = -(n*((-x(i))^(n-1)))./(a^n) ;
                fy(i) = -(n*((-y(i))^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*((-x(i))^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*((-y(i))^(n-2)))./(b^n) ;
                % kai_t =   atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                % kai_t =  2*pi -atan2((b^n*(-x(i)).^(n-1)),(a^n*(-y(i)).^(n-1)));               
                % kai_t(i) =  2*pi -atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   2*pi -atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            elseif (x(i)>=0 && y(i)<=0)
                beta(i) = x(i).^n /a^n + (-y(i)).^n /b^n ;
                fx(i) = (n*(x(i)^(n-1)))./(a^n) ;
                fy(i) = -(n*((-y(i))^(n-1)))./(b^n) ;
                fxx(i) = (n*(n-1)*(x(i)^(n-2)))./(a^n) ;
                fyy(i) = (n*(n-1)*((-y(i))^(n-2)))./(b^n) ;
                % kai_t =   -atan2(-(b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));                
                % kai_t =   atan2((b^n*(x(i)).^(n-1) ),(a^n*(-y(i)).^(n-1)));           
                % kai_t(i) =   atan2((b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   atan2(abs(fx(i)),abs(fy(i)));
                kait_dot(i) = -((fy(i).*fxx(i).*xdot(i) - fx(i).*fyy(i).*ydot(i))./fy(i)^2).*(cos(kai_t(i))).^2 ;
                beta_dot(i) = fx(i).*xdot(i) + fy(i).*ydot(i) ;
                kaio_dot(i) = ((2*k_se)./((1 + k_se*(beta(i) - 1).^2).*sqrt((k_se*(beta(i) - 1)).^2 + 2*k_se))).*beta_dot(i);
            end

            kai_o(i) = pi/2 -  asin(1./(1 + k_se*(beta(i)-1).^2))  ;

            
            if beta(i) <=1
                kaid(i) = kai_t(i)  -  kai_o(i)  ;
                kaid_dot(i) = kait_dot(i) - kaio_dot(i) ;

            else
                 kaid(i) = kai_t(i)  +  kai_o(i)  ;
                 kaid_dot(i) = kait_dot(i) + kaio_dot(i) ;

            end

        end
%     end
end
end

function kai = fun_getkai(x,y)
global  a b n k_se clockwise xc yc

%------------------------ initial kai  ------------------------------%
for i = 1: length(x)
%     for j = 1:length(X)
        %--------------------   vector field  ----------------------------%
        if clockwise ==0
         if (x(i)>=0 && y(i)>=0)
                beta = x(i).^n /a^n + y(i).^n /b^n ;
%                 des_angle(i,j) = (pi + atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
                 kai_t(i) = (pi - atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
            elseif (x(i)<=0 && y(i)>=0)
                beta = (-x(i)).^n /a^n + y(i).^n /b^n ;
%                 kai_t =  -pi - atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                  kai_t(i) =  pi + atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
%                   kai_t =  pi + atan2((b^n*(-x(i)).^(n-1)),(a^n*(y(i)).^(n-1)));
            elseif (x(i)<=0 && y(i)<=0)
                beta = (-x(i)).^n /a^n + (-y(i)).^n /b^n ;
%                 kai_t =   atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                  kai_t(i) =  2*pi -atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
%                  kai_t =  2*pi -atan2((b^n*(-x(i)).^(n-1)),(a^n*(-y(i)).^(n-1)));
            elseif (x(i)>=0 && y(i)<=0)
                beta = x(i).^n /a^n + (-y(i)).^n /b^n ;
%                 kai_t =   -atan2(-(b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   atan2((b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
%                 kai_t =   atan2((b^n*(x(i)).^(n-1) ),(a^n*(-y(i)).^(n-1)));
            end
            kai_o(i) = pi/2 -  asin(1./(1 + k_se*(beta-1).^2))  ;



            if beta <=1
                kai(i) = kai_t(i)  -  kai_o(i)  ;
%                 kaid_dot(i) = kait_dot(i) - kaio_dot(i) ;


            else
                 kai(i) = kai_t(i)  +  kai_o(i)  ;
%                  kaid_dot(i) = kait_dot(i) + kaio_dot(i) ;


            end

        else
                      if (x(i)>=0 && y(i)>=0)
                beta(i) = x(i).^n /a^n + y(i).^n /b^n ;
%                 kai_t = (pi + atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
                 kai_t(i) = (pi - atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1))));
            elseif (x(i)<=0 && y(i)>=0)
                beta(i) = (-x(i)).^n /a^n + y(i).^n /b^n ;
%                 kai_t =  -pi - atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                  kai_t(i) =  pi + atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
%                   kai_t =  pi + atan2((b^n*(-x(i)).^(n-1)),(a^n*(y(i)).^(n-1)));
            elseif (x(i)<=0 && y(i)<=0)
                beta(i) = (-x(i)).^n /a^n + (-y(i)).^n /b^n ;
%                 kai_t =   atan2(-(b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
                  kai_t(i) =  2*pi -atan2((b^n*abs(x(i)).^(n-1)),(a^n*abs(y(i)).^(n-1)));
%                  kai_t =  2*pi -atan2((b^n*(-x(i)).^(n-1)),(a^n*(-y(i)).^(n-1)));
            elseif (x(i)>=0 && y(i)<=0)
                beta(i) = x(i).^n /a^n + (-y(i)).^n /b^n ;
%                 kai_t =   -atan2(-(b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
                kai_t(i) =   atan2((b^n*abs(x(i)).^(n-1) ),(a^n*abs(y(i)).^(n-1)));
%                 kai_t =   atan2((b^n*(x(i)).^(n-1) ),(a^n*(-y(i)).^(n-1)));
            end
            kai_o(i) = pi/2 -  asin(1./(1 + k_se*(beta(i)-1).^2))  ;


            if beta(i) <=1
                kai(i) = kai_t(i)  -  kai_o(i)  ;
%                 kaid_dot(i) = kait_dot(i) - kaio_dot(i) ;

            else
                 kai(i) = kai_t(i)  +  kai_o(i)  ;
%                  kaid_dot(i) = kait_dot(i) + kaio_dot(i) ;
            end

        end
%     end
end
end


function [beta, fx, fxx,fy,fyy] = get_beta(x,y)
global a b n xc yc
for i = 1:length(x)
    if ((x(i,1)-xc)>=0) && ((y(i,1)-yc)>=0)
        beta = ((x(i,1) - xc).^n)./a^n + ((y(i,1)-yc).^n)./b^n;
        fx = (n*(x(i,1) - xc).^(n-1))./a^n ;
        fxx = (n*(n-1)*(x(i,1) - xc).^(n-2))./a^n ;
        fy = (n*(y (i,1)- yc).^(n-1))./b^n ;
        fyy = (n*(n-1)*(y(i,1) - yc).^(n-2))./b^n;
    elseif ((x(i,1)-xc)<0) && ((y(i,1)-yc)>0)
        beta = (-(x(i,1) - xc).^n)./a^n + ((y(i,1)-yc).^n)./b^n;
        fx = -(n*(x(i,1) - xc).^(n-1))./a^n ;
        fxx = -(n*(n-1)*(x(i,1) - xc).^(n-2))./a^n ;
        fy = (n*(y(i,1) - yc).^(n-1))./b^n ;
        fyy = (n*(n-1)*(y(i,1) - yc).^(n-2))./b^n;
    elseif ((x(i,1)-xc)<0) &&  ((y(i,1)-yc)<0)
        beta = ((x(i,1) - xc).^n)./a^n - ((y(i,1)-yc).^n)./b^n;
        fx = -(n*(x(i,1) - xc).^(n-1))./a^n ;
        fxx = -(n*(n-1)*(x(i,1) - xc).^(n-2))./a^n ;
        fy = -(n*(y(i,1) - yc).^(n-1))./b^n ;
        fyy = -(n*(n-1)*(y(i,1) - yc).^(n-2))./b^n;
    elseif ((x(i,1)-xc)>0) &&  ((y(i,1)-yc)<0)
        beta = ((x(i,1) - xc).^n)./a^n - ((y(i,1)-yc).^n)./b^n;
        fx = (n*(x(i,1) - xc).^(n-1))./a^n ;
        fxx = (n*(n-1)*(x(i,1) - xc).^(n-2))./a^n ;
        fy = -(n*(y(i,1) - yc).^(n-1))./b^n ;
        fyy = -(n*(n-1)*(y(i,1) - yc).^(n-2))./b^n;
    end
end
end

function path_error = fun_error(x,y)
global xc yc a b n
for i = 1:length(x)
    beta(i,1) = abs((x(i,1) - xc).^n)./a^n + abs((y(i,1)-yc).^n)./b^n;
%     [beta(i,1), fx(i,1), fxx(i,1),fy(i,1),fyy(i,1)] = get_beta(x(i,1),y(i,1));
    path_error(i,1) = beta(i,1) - 1;
end
end

function  [kaid0_ell,kaidot0_ell]= get_intialellipse(x0,y0)
global a_ell b_ell Vg k_se

beta_ell0 = (x0^2)./a_ell^2 + (y0^2)./b_ell^2 ;
fx_e0 = (2*x0)./a_ell^2 ;
fxx_e0 = 2./a_ell^2 ;
fy_e0 = (2*y0)./b_ell^2 ;
fyy_e0 = 2./b_ell^2 ;
kait_e0 =  atan2((b_ell^2*x0 ),(-a_ell^2*y0));

kaio_e0 = pi/2 -  asin(1./(1 + k_se*(beta_ell0-1).^2))  ;
if beta_ell0 <=1
    kaid0_ell = kait_e0 - kaio_e0 ;
    xdot0 = Vg*cos(kaid0_ell);
    ydot0  = Vg*sin(kaid0_ell);
    kaitdot_ell0 = -((fy_e0*fxx_e0*xdot0  - fx_e0*fyy_e0*ydot0 )./fy_e0^2)*(cos(kait_e0))^2 ;
    betadot_ell0 = fx_e0*xdot0  + fy_e0*ydot0  ;
    kaiodot_ell0 = ((2*k_se)./((1 + k_se*(beta_ell0 - 1)^2)*sqrt((k_se*(beta_ell0 - 1))^2 + 2*k_se)))*betadot_ell0;
    kaidot0_ell = kaitdot_ell0 - kaiodot_ell0;
else
    kaid0_ell = kait_e0 + kaio_e0 ;
    xdot0  = Vg*cos(kaid0_ell);
    ydot0  = Vg*sin(kaid0_ell);
    kaitdot_ell0 = -((fy_e0*fxx_e0*xdot0  - fx_e0*fyy_e0*ydot0)./fy_e0^2)*(cos(kait_e0))^2 ;
    betadot_ell0 = fx_e0*xdot0  + fy_e0*ydot0  ;
    kaiodot_ell0 = ((2*k_se)./((1 + k_se*(beta_ell0 - 1)^2)*sqrt((k_se*(beta_ell0 - 1))^2 + 2*k_se)))*betadot_ell0;
    kaidot0_ell = kaitdot_ell0 + kaiodot_ell0;
end
end

function [kaid0,kaid_dot0] = get_initialvalue(x0,y0)
global a b n Vg k_se

if (x0>=0 && y0>=0)
                beta0 = x0.^n /a^n + y0.^n /b^n ;
                fx0 = (n*(x0^(n-1)))./(a^n) ;
                fy0 = (n*(y0^(n-1)))./(b^n) ;
                fxx0 = (n*(n-1)*(x0^(n-2)))./(a^n) ;
                fyy0 = (n*(n-1)*(y0^(n-2)))./(b^n) ;
                % kai_t = (pi + atan2(-(b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1))));                
                % kait0 = (pi - atan2((b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1))));
                kait0 =   pi -atan2(abs(fx0),abs(fy0));
                
            elseif (x0<=0 && y0>=0)
                beta0 = (-x0).^n /a^n + y0.^n /b^n ;
                fx0 = -(n*((-x0)^(n-1)))./(a^n) ;
                fy0 = (n*(y0^(n-1)))./(b^n) ;
                fxx0 = (n*(n-1)*((-x0)^(n-2)))./(a^n) ;
                fyy0 = (n*(n-1)*(y0^(n-2)))./(b^n) ;
                % kai_t =  -pi - atan2(-(b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1))); 
                % kai_t =  pi + atan2((b^n*(-x0).^(n-1)),(a^n*(y0).^(n-1)));               
                % kait0 =  pi + atan2((b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1)));
                kait0 =   pi + atan2(abs(fx0),abs(fy0));
%                 kait_dot0 = -((fy0.*fxx0.*xdot0 - fx0.*fyy0.*ydot0)./fy0^2).*(cos(kait0)).^2 ;
%                 beta_dot0 = fx0.*xdot0 + fy0.*ydot0 ;
%                 kaio_dot0 = ((2*k_se)./((1 + k_se*(beta0 - 1).^2).*sqrt((k_se*(beta0 - 1)).^2 + 2*k_se))).*beta_dot0;
            elseif (x0<=0 && y0<=0)
                beta0 = (-x0).^n /a^n + (-y0).^n /b^n ;
                fx0 = -(n*((-x0)^(n-1)))./(a^n) ;
                fy0 = -(n*((-y0)^(n-1)))./(b^n) ;
                fxx0 = (n*(n-1)*((-x0)^(n-2)))./(a^n) ;
                fyy0 = (n*(n-1)*((-y0)^(n-2)))./(b^n) ;
                % kai_t =   atan2(-(b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1)));
                % kai_t =  2*pi -atan2((b^n*(-x0).^(n-1)),(a^n*(-y0).^(n-1)));               
                % kait0 =  2*pi -atan2((b^n*abs(x0).^(n-1)),(a^n*abs(y0).^(n-1)));
                kait0 =   2*pi -atan2(abs(fx0),abs(fy0));
%                 kait_dot0 = -((fy0.*fxx0.*xdot0 - fx0.*fyy0.*ydot0)./fy0^2).*(cos(kait0)).^2 ;
%                 beta_dot0 = fx0.*xdot0 + fy0.*ydot0 ;
%                 kaio_dot0 = ((2*k_se)./((1 + k_se*(beta0 - 1).^2).*sqrt((k_se*(beta0 - 1)).^2 + 2*k_se))).*beta_dot0;
            elseif (x0>=0 && y0<=0)
                beta0 = x0.^n /a^n + (-y0).^n /b^n ;
                fx0 = (n*(x0^(n-1)))./(a^n) ;
                fy0 = -(n*((-y0)^(n-1)))./(b^n) ;
                fxx0 = (n*(n-1)*(x0^(n-2)))./(a^n) ;
                fyy0 = (n*(n-1)*((-y0)^(n-2)))./(b^n) ;
                % kai_t =   -atan2(-(b^n*abs(x0).^(n-1) ),(a^n*abs(y0).^(n-1)));                
                % kai_t =   atan2((b^n*(x0).^(n-1) ),(a^n*(-y0).^(n-1)));           
                % kait0 =   atan2((b^n*abs(x0).^(n-1) ),(a^n*abs(y0).^(n-1)));
                kait0 =   atan2(abs(fx0),abs(fy0));
%                 kait_dot0 = -((fy0.*fxx0.*xdot0 - fx0.*fyy0.*ydot0)./fy0^2).*(cos(kait0)).^2 ;
%                 beta_dot0 = fx0.*xdot0 + fy0.*ydot0 ;
%                 kaio_dot0 = ((2*k_se)./((1 + k_se*(beta0 - 1).^2).*sqrt((k_se*(beta0 - 1)).^2 + 2*k_se))).*beta_dot0;
            end

            kaio0 = pi/2 -  asin(1./(1 + k_se*(beta0-1).^2))  ;

            
            if beta0 <=1
                kaid0 = kait0  -  kaio0  ;
                xdot0 = Vg*cos(kaid0);
                ydot0 = Vg*sin(kaid0);
                kait_dot0 = -((fy0.*fxx0.*xdot0 - fx0.*fyy0.*ydot0)./fy0^2).*(cos(kait0)).^2 ;
                beta_dot0 = fx0.*xdot0 + fy0.*ydot0 ;
                kaio_dot0 = ((2*k_se)./((1 + k_se*(beta0 - 1).^2).*sqrt((k_se*(beta0 - 1)).^2 + 2*k_se))).*beta_dot0;
                kaid_dot0 = kait_dot0 - kaio_dot0 ;

            else
                 kaid0 = kait0  +  kaio0  ;
                 xdot0 = Vg*cos(kaid0);
                 ydot0 = Vg*sin(kaid0);
                 kait_dot0 = -((fy0.*fxx0.*xdot0 - fx0.*fyy0.*ydot0)./fy0^2).*(cos(kait0)).^2 ;
                 beta_dot0 = fx0.*xdot0 + fy0.*ydot0 ;
                 kaio_dot0 = ((2*k_se)./((1 + k_se*(beta0 - 1).^2).*sqrt((k_se*(beta0 - 1)).^2 + 2*k_se))).*beta_dot0;
                 kaid_dot0 = kait_dot0 + kaio_dot0 ;

            end
end
