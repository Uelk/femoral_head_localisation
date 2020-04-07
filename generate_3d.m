%% Visualize 3D Data

% profile on;  %profile viewer   %Optimize Code
close all; clc; clear;

%% Settings: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bool_spehre             = 1;      % Make Sphere Visible?
Overlay_multiple_data   = 1;      % How many data sets?
bool_plot_point_cloud   = 1;      % Plot every individual w_t?
bool_ratio_spehre       = 1;      % Show the Threshold of 2.5mm for pointcloud
trusted_radius          = 2.5e-3; % Trusted radius for deviation from Barycenter in meters

bool_import_theta_u     = 0;      % Set whether the imported data has to be converted from theta_u to rotation matrix first

% File import settings
FolderName              = 'C:\Users\Phrittich\Desktop\Praktikum Aesculap\3D Visualization\data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Einheitsvektor
I = [1 0 0  ;
    0 1 0  ;
    0 0 1 ];


for struct_index = 1 : Overlay_multiple_data
    [FileName, FolderName]  = uigetfile(fullfile( FolderName,'*.txt'));
    FullPath                = fullfile(FolderName, FileName);
    if(bool_import_theta_u ~= 1)
        imported_data                      = func_import_xyz_rot_mat(FullPath, [4, Inf]);  %Format: X Y Z R11, R12,R13, R21,....
        struct_data(struct_index).data     = imported_data;
        struct_data(struct_index).filename = FileName;
        
    else
        imported_data                      = func_import_xyz_theta_u(FullPath, [4, Inf]);  %Format: X Y Z Theta11 Theta12 Theta 13
        struct_data(struct_index).filename = FileName;
        
        % Conversion from theta_u to rotation matrix:
        
        % Rotation Representation
        % axis/Angle : (u_vec , theta)
        % Theta_u    : theta * u_vec
        % Matrix     : [r1 r2 r3 ; ....]
        
        theta_u = imported_data(:,4:6);
        for i=1 : (size(theta_u,1))                             
            theta = norm(theta_u(i,1:3));                       % angle         
            u_vec = theta_u(i,1:3) / norm(theta_u(i,1:3));      % axis
            
            cross_u_vec =  [    0    , -u_vec(3) ,   u_vec(2);
                             u_vec(3),     0     ,  -u_vec(1);
                            -u_vec(2), u_vec(1)  ,      0   ];
            
            cross_u_vec2 = cross_u_vec * cross_u_vec ;
                        
            % Rodrigues-Formel: Matrix exponential of rotation
            Rotation_matrix = I + sin(theta)* cross_u_vec + (1-cos(theta)) * cross_u_vec2 ;

            struct_data(struct_index).data(i,4:12)    = [Rotation_matrix(1,1) Rotation_matrix(1,2) Rotation_matrix(1,3) Rotation_matrix(2,1) Rotation_matrix(2,2) Rotation_matrix(2,3) Rotation_matrix(3,1) Rotation_matrix(3,2) Rotation_matrix(3,3)];
        end
        struct_data(struct_index).data(:,1:3)     = imported_data(:,1:3);
    end
end




%% AOS - Algebraic One Step

for calc_index = 1 : Overlay_multiple_data
    disp(sprintf('\n\n ########## Results for File [%s]: ################\n' ,struct_data(calc_index).filename));            % Show progress during runtime
    R = zeros(size( struct_data(calc_index).data, 1) , 6); % Initialize
    t = zeros(size( struct_data(calc_index).data, 1) , 1); % Initialize
    R_index = 1;
    t_index = 1;
    
    % Set up the R_i and t_i matrices and already place them beneath eachother as seen in AOS-Formular
    for data_point =  1 :  size(struct_data(calc_index).data,1)
        R(R_index, 1:6)         = [struct_data(calc_index).data(data_point,4:6)     -I(1,1:3)];
        R_index                 = R_index + 1;
        R(R_index, 1:6)         = [struct_data(calc_index).data(data_point,7:9)     -I(2,1:3)];
        R_index                 = R_index + 1;
        R(R_index, 1:6)         = [struct_data(calc_index).data(data_point,10:12)   -I(3,1:3)];
        R_index                 = R_index + 1;
        
        t(t_index)              = struct_data(calc_index).data(data_point,1);
        t_index                 = t_index + 1;
        t(t_index)              = struct_data(calc_index).data(data_point,2);
        t_index                 = t_index + 1;
        t(t_index)              = struct_data(calc_index).data(data_point,3);
        t_index                 = t_index + 1;
    end
    
    %      A * x = b
    %     (A^T * A) * x = A^T * b
    %      x = (A^T * A)^-1 * A^T * b
    A                                   = R;
    AA                                  = A' * A;
    result                              = inv(AA) * (A' * -t);
    
    DRF_t                               = result(1:3);
    struct_data(calc_index).DRF         = DRF_t;
    
    w_t                                 = result(4:6);
    struct_data(calc_index).w_t_algorith= w_t;
    
    fprintf('w_t from the AOS-Algorithm = [%0.5d  %0.5d  %0.5d] \n\n' , w_t(1), w_t(2),w_t(3) ) ;
    
    
    % Calc every individual w_t from the calculated DFT_t to compare to w_t_algorith :  w_t = R_i * DRF_t + t_i
    
    struct_index_individ = 1;
    for individ_index =  1 : 3 : size(R)
        R_i = [ R(individ_index, 1:3) ; R(individ_index+1, 1:3) ; R(individ_index+2, 1:3) ];
        t_i = [ t(individ_index, 1)   ; t(individ_index+1, 1)   ; t(individ_index+2, 1)   ];
        w_t_individ = R_i * DRF_t + t_i;
        struct_data(calc_index).w_t_individ(:,struct_index_individ) = w_t_individ;
        struct_index_individ = struct_index_individ + 1;
    end
    
    
    
    
    %Calculate the barycenter, which is the mean of all the individual w_t
    struct_data(calc_index).barycenter(1,1) = mean( struct_data(calc_index).w_t_individ(1,:) );
    struct_data(calc_index).barycenter(2,1) = mean( struct_data(calc_index).w_t_individ(2,:) );
    struct_data(calc_index).barycenter(3,1) = mean( struct_data(calc_index).w_t_individ(3,:) );
    fprintf('Barycenter from the mean of all individual w_t = [%0.5d  %0.5d  %0.5d] \n\n' , struct_data(calc_index).barycenter(1),struct_data(calc_index).barycenter(2), struct_data(calc_index).barycenter(3) ) ;
    
    
    
    % Difference between Barycenter and w_t_algorith
    
    diff = w_t - struct_data(calc_index).barycenter;
    diff_norm = norm(diff);
    fprintf('Difference of w_t to Barycenter = [%0.5d  %0.5d  %0.5d] \n\n', diff(1), diff(2),diff(3) ) ;
    fprintf('Norm of Difference = %0.5d \n\n' ,diff_norm ) ;
    
    
    
    % Calculate RMSE of all the individual w_t to the Barycenter
    
    error           = (struct_data(calc_index).w_t_individ - struct_data(calc_index).barycenter);        % Errors
    erro_sqr        = (struct_data(calc_index).w_t_individ - struct_data(calc_index).barycenter).^2 ;    % Squared Error
    mean_sqr_error  = [ mean(erro_sqr(1,:)) ; mean(erro_sqr(2,:)) ; mean(erro_sqr(3,:)) ]  ;             % Mean Squared Error
    RMSE            = sqrt(mean_sqr_error);                                                              % Root Mean Squared Error
    
    fprintf('The RMSE of all the individual w_t to the Barycenter = [%0.5d  %0.5d  %0.5d] \n\n' ,RMSE(1) ,RMSE(2), RMSE(3) ) ;
    
    %
    %     msgbox( { sprintf('The Barycenter for File: %s was accepted', struct_data(calc_index).filename) , sprintf('The Barycenter = [%0.5d  %0.5d  %0.5d]', struct_data(calc_index).barycenter(1),struct_data(calc_index).barycenter(2), struct_data(calc_index).barycenter(3)) },'Accepted','help');
    %
    %     msgbox( { sprintf('The Barycenter for File: %s was NOT accepted', struct_data(calc_index).filename) , sprintf('The Barycenter = [%0.5d  %0.5d  %0.5d]', struct_data(calc_index).barycenter(1),struct_data(calc_index).barycenter(2), struct_data(calc_index).barycenter(3)) },'Not Accepted','error');
    %
    
    
end


%% Plotting
for plot_index = 1 : Overlay_multiple_data
    figure(1);
    
    if (plot_index == 1)        %Set different color for different datasets
        RGB = [256 0 0]/256 ;   % First Always Red!
    else
        r = rand(1,3); % Random colour
        RGB = [256*r(1) 256*r(2) 265*r(3)]/256 ;
        RGB = [0 0 256]/256 ;
    end
    
    
    % Scatter Plot for the measured points on the sphere
    %                          X                                    Y                                   Z                      Size Color  Form          Legend-Name
    scatt = scatter3(struct_data(plot_index).data(:,1) , struct_data(plot_index).data(:,2) , struct_data(plot_index).data(:,3), 20 , RGB , '.' ,'DisplayName', sprintf('3D-Data for File: %s',struct_data(plot_index).filename));
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on;
    axis equal;
    hold on;
    
    
    
    %Plot w_t -> diamond
    p = scatter3(struct_data(plot_index).w_t_algorith(1), struct_data(plot_index).w_t_algorith(2) , struct_data(plot_index).w_t_algorith(3) , 50 , RGB, 'd','DisplayName', sprintf('w_t from algorithm for file: %s',struct_data(plot_index).filename));  % d = diamond  % 50 = size
    text_handle = text(struct_data(plot_index).w_t_algorith(1), struct_data(plot_index).w_t_algorith(2) , struct_data(plot_index).w_t_algorith(3), ['\leftarrow' , sprintf('w_{t}_{mean} %d',plot_index)], 'Color',RGB,'FontSize',7);
    hold on;
    
    %Plot barycenter -> pentagram
    p = scatter3(struct_data(plot_index).barycenter(1), struct_data(plot_index).barycenter(2) , struct_data(plot_index).barycenter(3) , 50 , RGB, 'p','DisplayName', sprintf('Barycenter from mean of all w_t_individ for file: %s',struct_data(plot_index).filename));  % d = diamond  % 50 = size
    text_handle = text(struct_data(plot_index).barycenter(1), struct_data(plot_index).barycenter(2) , struct_data(plot_index).barycenter(3), ['\leftarrow' , sprintf('BC %d',plot_index)], 'Color',RGB,'FontSize',7);
    hold on;
    
    
    % Plot point cloud for each individual w_t
    if(bool_plot_point_cloud ==1)
        scatter3(struct_data(plot_index).w_t_individ(1,:), struct_data(plot_index).w_t_individ(2,:) , struct_data(plot_index).w_t_individ(3,:) , 20 , RGB, '.','HandleVisibility','off'); %HandleVis -> Off: No Legend!
    end
    
    
    %Sphere Plot of error Threshold 2.5mm
    if(bool_ratio_spehre == 1)
        figure(1);
        [X,Y,Z] = sphere(100); %Sphere with 100 Faces
        sphe = surf(X*trusted_radius + struct_data(plot_index).barycenter(1), Y*trusted_radius + struct_data(plot_index).barycenter(2), Z*trusted_radius + struct_data(plot_index).barycenter(3), 'FaceAlpha',0 , 'EdgeAlpha',0.07,'HandleVisibility','off');
        sphe.FaceColor = RGB;
        sphe.EdgeColor = RGB;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        grid on;
        axis equal;
    end
    
    
end

%Sphere Plot of 1st File
if(bool_spehre == 1)
    
    figure(1);
    [X,Y,Z] = sphere(100); %Sphere with 100 Faces
    sphe = surf(X*norm(DRF_t) + struct_data(calc_index).w_t_algorith(1), Y*norm(DRF_t) + struct_data(calc_index).w_t_algorith(2), Z*norm(DRF_t) + struct_data(calc_index).w_t_algorith(3), 'FaceAlpha',0, 'EdgeAlpha',0.07,'HandleVisibility','off');
    sphe.FaceColor = [1 0 0];
    axis equal;
else % Create a Sphere with no outlines
    figure(1);
    [X,Y,Z] = sphere(100); %Sphere with 100 Faces
    sphe = surf(X*norm(DRF_t) + struct_data(calc_index).w_t_algorith(1), Y*norm(DRF_t) + struct_data(calc_index).w_t_algorith(2), Z*norm(DRF_t) + struct_data(calc_index).w_t_algorith(3), 'FaceAlpha',0, 'EdgeAlpha',0,'HandleVisibility','off');
    sphe.FaceColor = [1 0 0];
    axis equal;
end


title_handle        = title('3D-Data for Files');
title_handle.Color  = [0.9290 0.6940 0.1250];
set(title_handle,'Interpreter','none')      %Damit die Unterstriche nicht den Titel verkacken
lgnd                = legend('Location','southeast');
set(lgnd,'Interpreter','none')              %Damit die Unterstriche nicht die Legende verkacken

