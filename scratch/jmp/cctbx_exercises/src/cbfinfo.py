"""A module to print information about a CBF file"""

import pycbf

def print_info(cbf_path):
    """Print out a load of data held in the CBF file. 
    
    This is by no means a full list of the data contained in the file, it's 
    mainly for debugging and development purposes. The data that will be 
    printed is the following:
    - The number of categories and the name of each category
    - The number of rows and columns and the name of each column
    - The type of each element in each row/column element.
    
    :param cbf_path: The path to the cbf file
    """

    # Read the CBF file    
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cbf_path, pycbf.MSG_DIGEST)
    cbf_handle.rewind_datablock()
    
    # Select the first datablock and rewind all the categories
    cbf_handle.select_datablock(0)
    cbf_handle.rewind_category()

    # Count the number of categories and loop through them
    num_categories = cbf_handle.count_categories()
    for i in range(num_categories):

        # Select the ith category and print its name
        cbf_handle.select_category(i)
        category_name = cbf_handle.category_name()
        print "Category:", i, category_name    

        # Count the number of rows and columns in the category
        # and print them
        num_rows = cbf_handle.count_rows()
        num_cols = cbf_handle.count_columns()
        print "\tNum (rows, cols)", (num_rows, num_cols)

        # Rewind the columns and print the name of each
        cbf_handle.rewind_column()
        for i in range(num_cols):
            cbf_handle.select_column(i)
            column_name = cbf_handle.column_name()
            print '\tColumn:', i, column_name

        # Loop through all rows and columns and print the
        # type of the data stored in that table element
        for j in range(num_rows):
            cbf_handle.select_row(j)
            cbf_handle.rewind_column()
            print '\t\tRow:', j, cbf_handle.get_value()
            for i in range(num_cols):
                cbf_handle.select_column(i)
                type_of_value = cbf_handle.get_typeofvalue()
                if type_of_value.find('dblq') > -1:
                    value = cbf_handle.get_value()
                elif type_of_value.find('text') > -1:
                    value = cbf_handle.get_value()
                    value = '\n\t\t\t'.join(value.split('\n'))
                elif type_of_value.find('word') > -1:
                    value = cbf_handle.get_value()
                elif type_of_value.find('sglq') > -1:
                    value = cbf_handle.get_value()
                else:
                    value = '...'
                print "\t\tColumn", i, "Type:", type_of_value, value
  

                