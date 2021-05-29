function savetxt(filename,x_axis,y_axis,matrix)
    size_ = size(matrix);
    ii = size_(1);
    jj = size_(2);
    outfile = fopen(filename,'w+');
    for i = 1:ii
        for j = 1:jj;
            fprintf(outfile,'%15.10f  %15.10f %20.15f',x_axis(i),y_axis(j),matrix(i,j));
            fprintf(outfile,'\n');
        end;
    end;
fclose(outfile);          
end
    