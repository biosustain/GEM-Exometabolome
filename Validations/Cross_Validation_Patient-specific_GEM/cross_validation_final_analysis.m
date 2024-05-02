table_1 = load("Cross_validation_table_1-5.mat");
table_2 = load("Cross_validation_table_6-10.mat");
table_3 = load("Cross_validation_table_11-15.mat");
table_4 = load("Cross_validation_table_16-20.mat");
table_5 = load("Cross_validation_table_21-25.mat");
table_6 = load("Cross_validation_table_26-30.mat");
table_7 = load("Cross_validation_table_31-35.mat");
table_8 = load("Cross_validation_table_36-40.mat");
table_9 = load("Cross_validation_table_41-45.mat");
table_10 = load("Cross_validation_table_46-50.mat");
table_11 = load("Cross_validation_table_51-55.mat");
table_12 = load("Cross_validation_table_56-60.mat");
table_13 = load("Cross_validation_table_61-65.mat");
table_14 = load("Cross_validation_table_66-70.mat");
table_15 = load("Cross_validation_table_71-75.mat");
table_16 = load("Cross_validation_table_76-80.mat");
table_17 = load("Cross_validation_table_81-85.mat");
table_18 = load("Cross_validation_table_86-90.mat");
table_19 = load("Cross_validation_table_91-95.mat");

% Initialize the final cell array
final_table = cell(53, 95);

% Loop over all the table structures
for i = 1:19
    % Get the table structure
    table = eval(['table_', num2str(i)]);
    
    % Get the cross_validation_table from the structure
    cross_validation_table = table.cross_validation_table;
    
    % Get the columns to be copied (5 columns for each table)
    cols = (i-1)*5 + 1 : i*5;
    
    % Copy the non-empty columns to the final table
    final_table(:, cols) = cross_validation_table(:, cols);
end


for i = 1:53
    for j = 1:95
        if abs(final_table{i,j}(5)) + abs(final_table{i,j}(5)) == 0
            final_table{i,j}(1) = 0;
        end
    end
end
