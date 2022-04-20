import pandas as pd
import openpyxl as pxl
import os

# ------------------------------------------------------

# Save the input dataframe to the specified sheet name of filename file
def saveAsNewSheetToExistingFile(filename,newDF,newSheetName):
    if os.path.exists(filename):

        excel_book = pxl.load_workbook(filename)

        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            # Your loaded workbook is set as the "base of work"
            writer.book = excel_book

            # Loop through the existing worksheets in the workbook and map each title to\
            # the corresponding worksheet (that is, a dictionary where the keys are the\
            # existing worksheets' names and the values are the actual worksheets)
            writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}

            # Write the new data to the file without overwriting what already exists
            newDF.to_excel(writer, newSheetName)

            # Save the file
            writer.save()
    else:
        newDF.to_excel(filename, newSheetName)
        
    return


# ------------------------------------------------------

# saveDF_to_CSV_GZ_no_timestamp
def saveDF_to_CSV_GZ_no_timestamp(df,filename):
    from gzip import GzipFile
    from io import TextIOWrapper
    with TextIOWrapper(GzipFile(filename, 'w', mtime=0), encoding='utf-8') as fd:
        df.to_csv(fd,index=False,compression='gzip')
        
    return