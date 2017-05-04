def remove_lines(in_file, out_file, match_string, above, below):
    """
    Read from input file, search for lines containing a string, remove matching lines and a number of lines above
    and below that line, and write out remaining lines to output file.

    :param in_file: open file object to read
    :param out_file: open file object to write
    :param match_string: str -- string to search for
    :param above: int -- number of lines to remove above the matched line
    :param below: int -- number of lines to remove below the matched line
    """
    assert(isinstance(above, int))
    assert (isinstance(below, int))
    buff = []
    line = in_file.readline()
    while line:
        if match_string in line:
             buff = []
             for _ in range(below):
                 in_file.readline()
        else:
            if len(buff) == above:
                out_file.write(buff[0])
                buff = buff[1:]
            buff.append(line)
        line = in_file.readline()
    out_file.write(''.join(buff))

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print('Usage: {} match above below'.format(sys.argv[0]))
        print('  Remove every instance of matching string from stdin and a number of lines above and below the match.')
        sys.exit(-1)
    remove_lines(sys.stdin, sys.stdout, sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))