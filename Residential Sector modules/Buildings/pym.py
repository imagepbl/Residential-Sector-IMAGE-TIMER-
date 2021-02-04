import numpy as np
import pandas as pd
import re as re
import os.path

from csv import Sniffer as Sniffer
from string import whitespace as whitespace

LINESEP = "\n"

# ............................................................................ #
def process_header(line):
    """

    """
    raw_data = ""
    # > the header line contains "[d1,d2,...,dn](t) = [t1", with n [di]
    #   dimensions, where the time dimension "(t)" is optional.
    # >> define [endmarker] of variable dependent on time-dependency
    has_time = True if "(t)" in line else False
    
    # >> parse the [di] dimension values
    # TODO: check if this works
    dim_match = re.search(r"\[(.*?)\]", line).group(1)
    if dim_match is not None:
        # dim_match = dim_match.group(1)
        dim_match = re.split(",", dim_match)
        dims = [int(d) for d in dim_match]
    else:
        dims = [1]
    # >> the data values length includes 1 year value if [has_time]
    data_length = np.prod(dims) + 1 if has_time else np.prod(dims)
    
    # > the following "[" and "t1" are only there in the case of a
    #   time-dependent variable; t1 can also be on the next line(s).
    # >> parse the remainder of the line
    if has_time:
        start_array = line.rfind("[") + 1
        if line[start_array:].isdigit():
            raw_data = line[start_array:]

    return has_time, dims, data_length, raw_data
# ............................................................................ #


def read_mym(filename, path=""):
    """
    Read a MyM data file. Return numpy array.
    
    The MyM data file should contain a single variable, that can be
    time-dependent. 
    
    Parameters
    ----------
    filename : string
        name of the file to be read in
    path : string
        path to the filename, either relative or absolute
    
    Returns
    -------
    header : string
        contains comment lines and variable specification, the last line should
        end with '[a1,...,an](t) = [' or '[a1,...,an] ='
    data : np.array
        contains the data as a [t,a1,...,an] shaped numpy matrix
    time : np.array
        vector with time values, with time.shape = [t,1]
    """

    # Open file and read data
    filename = os.path.join(path, filename)
    
    with open(filename, "r") as mym_file:
        header = []
        line = mym_file.readline().strip(whitespace + ",")
        while line[0] == "!":
            header.append(line)
            line = mym_file.readline().strip(whitespace + ",")
        # > process header and put data on header line into string 
        has_time, dimensions, data_length, raw_data = process_header(line)
        content = mym_file.read().replace("\n", " ").rstrip("];" + whitespace)
    
    # Process data
    # > transform to numpy array, using the csv.Sniffer to find the delimiter.
    # > `first_chunk` is a large enough string for the sniffer to find
    #   the delimiter and small enough to be fast.
    first_chunk = 64
    delimiters = "".join([",", whitespace])
    mym_format = Sniffer().sniff(content[:first_chunk], delimiters=delimiters)
    raw_data = mym_format.delimiter.join([s for s in [raw_data, content] if s])
    raw_data = np.fromstring(raw_data, sep=mym_format.delimiter)

    # > find the desired dimensions, where `time_length` == 1 for
    #   time-independent data 
    if raw_data.size % data_length == 0:
        time_length = int(raw_data.size / data_length)
        raw_dimensions = (time_length,data_length)
        target_dimensions = tuple([time_length] + dimensions)
    else:
        raise RuntimeError("file dimensions are parsed incorrectly")

    # > reshape data in two steps to split off time vector, while reflecting
    #   original dimensions, where `time = data[:,0]`
    data = np.reshape(raw_data, raw_dimensions)
    if has_time:
        time, data = np.split(data,[1], axis=1)
        data = np.reshape(data, target_dimensions)
        time = np.squeeze(time)
        return data, time
    else:
        return data
# ............................................................................ #


def stringify(array, indent=4):
    """
    Generate data string from array. Return string.
    """
    formatter = lambda number: "{:14s}".format("{:.4f},".format(number))
    s = ""
    for row in array:
        line = [" " * indent]
        line += [formatter(x) for x in row]
        line += LINESEP
        s += "".join(line)
    
    return s
# ............................................................................ #


def get_yearly_data_table(data, dimensions, year, years, table):
    """
    Generate `data` for `year` in desired `table` format. Return numpy array.
    """
    if isinstance(data, pd.DataFrame):
        data_table = data.loc[:,[year]]
    elif isinstance(data, np.ndarray):
        data_table = data[year == years,...]
    
    data_table = shape_data_table(data_table, dimensions, table)
    
    return data_table
# ............................................................................ #


def shape_data_table(data_table, dimensions, table):
    """
    Shape `data_table` to desired table format. Return numpy array.
    """
    if isinstance(data_table, (pd.DataFrame, pd.Series)):
        data_table = data_table.values

    to_shape = {"wide": [-1, dimensions[-1]], "long": [-1, 1]}
    data_table = data_table.reshape(to_shape[table])

    return data_table


def print_mym(data, years=None, name="", table="long", comment=""):
    """
    Print `data` in MyM data format string. Return string.

    The time dimension in `years` is treated as a separate, by definition
    slowest changing dimension. Of the other `dimensions`, the leftmost
    dimension is the slowest changing dimension after `years`. The rightmost
    dimension is the fastest changing dimension.
    Long table format prints `data` in a single column, while wide table format
    prints `data` in matrix format with the fastest changing dimension on a
    single line.

    Parameters
    ----------
    data : pandas.DataFrame, pandas.Series or numpy.array
        a DataFrame is used for a time-dependent variable, where each column
        contains the data for a single year, for a time-independent variable a
        Series is expected. A numpy array can contain both types of variables.
    years : array_like
        contains time values associated with `data` when data is a numpy array,
        optional
    name : string
        name of the variable contained in `data`
    table : {"long", "wide"}
        ways of formatting MyM `data` in string
    comment : string
        description of `data`, optional

    Returns
    -------
    data_string : string
        `data` in MyM string format

    See also
    --------
    write_mym : Write MyM `data_string` to file
    read_mym : Read MyM data from file
    """
    if table not in ["long", "wide"]:
        table = "long"
        print("The [table] argument is not specified correctly, it should be"
              " either 'long' (default) or 'wide', using default.")

    # Generate [dimensions] and [datatype]
    # > [dimensions] is cast into a list for the correct string representation
    if isinstance(data, (pd.DataFrame, pd.Series)):
        try:
            # Clean [data] for known MultiIndex issues
            data.index = data.index.remove_unused_levels()
            if all([level in data.index.names for level in ["image_region","image_region_index"]]):
                data.index = data.index.droplevel(level="image_region_index")
            dimensions = list(data.index.levshape)
            #  Check dimensionality
            # > [data] should have a size that is equal to the product of its
            #   dimensions, otherwise [data] has duplicates or is incomplete.
            if not np.prod(dimensions) == data.shape[0]:
                raise Exception("Data for [{}] variable has an inconsistent"
                                " format. Expected [{}] index entries, got"
                                " [{}].".format(name, np.prod(dimensions), data.shape[0]))
        except:
            dimensions = [len(data.index)]
        datatype = data.values.dtype
        years = [year for year in data] if isinstance(data, pd.DataFrame) else None
    elif isinstance(data, np.ndarray):
        dimensions = list(data.shape[1:])
        datatype = data.dtype
    else:
        raise TypeError("Data for [{}] should be a numpy array or a pandas"
                        " DataFrame or Series, instead data is of {}.".format(name, type(data)))
    
    if not np.issubdtype(datatype, np.number):
        raise ValueError("Data for [{}] contains non-numeric entries.".format(name))

    # Print header
    if comment:
        comment = "!" + comment if comment[0] != "!" else comment
        comment = comment + LINESEP if comment[-1] != LINESEP else comment
    mym_datatype = "real" if np.issubdtype(datatype, np.floating) else "integer"
    time = "(t) = [" if years is not None else " ="
    dimension_string = str(dimensions) if np.prod(dimensions) > 1 else ""
    header = comment + "{} {}{}{}{}".format(mym_datatype, name, dimensions, time, LINESEP)

    # Print data
    data_string = [header]
    if years is None:
        data_string.append(stringify(shape_data_table(data, dimensions, table)))
    else:
        for year in years:
            data_string.append("{},{}".format(year, LINESEP))
            data_block = get_yearly_data_table(data, dimensions, year, years, table)
            data_string.append(stringify(data_block))

    # Close MyM data array
    # > last number in array should not be followed by a comma
    data_string[-1] = data_string[-1].rstrip("\n, ")
    if years is not None:
        data_string.append("]")
    data_string.append(";")

    return "".join(data_string)
# ............................................................................ #


def write_mym(data, years=None, table="long", variable_name="data",
        filename=None, filepath=None, comment=""):
    """
    Write mym `data_string` to file.
    """
    data_string = print_mym(data, years=years, name=variable_name, table=table,
                            comment=comment)

    filename = variable_name + ".dat" if not filename else filename
    filepath = os.path.join(filepath, filename) if filepath else filename
    with open(filename, "w+") as mym_file:
        mym_file.write(data_string)


if __name__ == "__main__":

    read_mym("floorspace.out")

    print("Performing pym integration tests.")
    dimensions = (4,5,6)
    data = np.random.rand(*dimensions)
    time = np.array([2015, 2020, 2025, 2030])
    variable_name = "test"
    filename = variable_name + ".dat"
    comment = ("This is a test data set of random numbers with dimensions {},"
               " where the first number indicates the number of time entries"
               .format(dimensions))

    print("  writing test data using [write_mym]:")
    try:
        write_mym(data, years=time, table="wide", variable_name=variable_name, 
                  comment=comment)
    except Exception as error:
        print("  x test data could not be written to mym-file, the following"
              " error occured:")
        print(error)
        raise error
    else:
        print("  > successfully written data to file")
    
    print("  reading test data using [read_mym]:")
    try:
        data_output, time_output = read_mym(filename)
    except Exception as error:
        print("Test data could not be read from [{}]".format(filename))
        print(error)
        raise error
    else:
        print("  > successfully read data into array")
    
    print("  checking equality of test input and output:")
    try:
        assert np.allclose(time_output, time, atol=1e-04)
        assert np.allclose(data_output, data, atol=1e-04)
    except AssertionError as error:
        print("  x time and data input and output do not match.")
        print(error)
        raise error
    else:
        print("  > input and output data match to required precision")
        print("Finished tests succesfully.")

    print()