function [] = part_2_b()
part_2_mesh_vs_current();
part_2_gap_vs_current();
part_2_sigma_vs_current();
end

function [] = part_2_mesh_vs_current()
    clear all;
    close all;
    clc;

    nx = 10:5:200;
    ny = nx;
    delta = 1;

    a_cond = 1e-2;
    b_cond = 1;

    L = nx*delta;
    W = ny*delta;

    box_w = 50;
    box_gap = 50;

    for i = 1:1:length(nx)
        [~,~,~,~,~,current_end(i),current_start(i)] = part_2_solver(nx(i),ny(i),delta/(nx(i)/nx(1)),a_cond,b_cond,box_w,box_gap);
    end
    
    figure(1);
    plot(nx,current_end,'b.-');hold on;
    plot(nx,current_start,'r.-');
    legend('Current at end of contact','Current at start of contact');
    xlabel('Mesh Size');
    ylabel('Current');    
    title('Assignment 2 Part 2 B - Current vs. Mesh Size with \Delta Adjustment');
    print(figure(1),'current_vs_mesh_adjusted_delta','-depsc');
end

function [] = part_2_gap_vs_current()
    clear all;
    close all;
    clc;
    
    nx = 100;
    ny = nx;
    delta = 1;

    a_cond = 1e-2;
    b_cond = 1;

    L = nx*delta;
    W = ny*delta;

    box_w = 20;
    box_gap = 1:50;

    for i = 1:1:length(box_gap)
        [~,~,~,~,~,current_end(i),current_start(i)] = part_2_solver(nx,ny,delta,a_cond,b_cond,box_w,box_gap(i));
    end
    
    figure(2);
    plot(box_gap,current_end,'b.-');hold on;
    plot(box_gap,current_start,'r.-');
    legend('Current at end of contact','Current at start of contact');
    xlabel('Gap Size');
    ylabel('Current');    
    title('Assignment 2 Part 2 B - Current vs. Gap Size with \Delta Adjustment');
    print(figure(2),'current_vs_box_width','-depsc');
end

function [] = part_2_sigma_vs_current()
    clear all;
    close all;
    clc;
    
    nx = 100;
    ny = nx;
    delta = 1;

    a_cond = 0:1e-2:1;
    b_cond = 1;

    L = nx*delta;
    W = ny*delta;

    box_w = 20;
    box_gap = 30;

    for i = 1:1:length(a_cond)
        [~,~,~,~,~,current_end(i),current_start(i)] = part_2_solver(nx,ny,delta,a_cond(i),b_cond,box_w,box_gap);
    end
    
    figure(3);
    plot(a_cond,current_end,'b.-');hold on;
    plot(a_cond,current_start,'r.-');
    legend('Current at end of contact','Current at start of contact');
    xlabel('\sigma(S)');
    ylabel('Current');    
    title('Assignment 2 Part 2 B - Current vs. \sigma with \Delta Adjustment');
    print(figure(3),'current_vs_sigma','-depsc');
end