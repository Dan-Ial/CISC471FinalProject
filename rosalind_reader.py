'''
    Simple file reader for taking in rosalind data, which will always
    have characters in a txt file as input separated by whitespace.
'''

def read_file(file_name):
    print("Reading file...")
    file_stream = open(file_name, "r")

    # read file data into data array
    try:
        data = file_stream.read().strip().split()
    except:
        print("File read failed. Please check that the file exists.")
        return None
    finally:
        # close file stream
        file_stream.close()
        print("File read successfully")
        return data