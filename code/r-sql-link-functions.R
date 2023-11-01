# Functions for establishing connections to databases and editing sql code
# Some based off of functions from gfdata in Utils.R
# db_connection() works with both 32bit and 64bit R


sf_db_connection <- function(server = "WDC-SQL2016-P\\SIOSP01") {
  require(DBI)
  DBI::dbConnect(odbc::odbc(),
                 driver = "SQL Server",
                 server = server)
}

mdb_connection <- function(db_file_path)  {
  require(DBI)
  # Make sure that the file exists before attempting to connect
  if (!file.exists(db_file_path)) {
    stop("DB file does not exist at ", db_file_path)
  }
  # Connect
  DBI::dbConnect(odbc::odbc(),
                 .connection_string = 
                   paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};", 
                          paste0("DBQ=", db_file_path)))
}

cap_words <- function( s, split="_" ) {
  cap <- function(i)
    paste(toupper(substring(i, 1, 1)),
          {i <- substring(i, 2); tolower(i)},
          sep = "", collapse = "_")
  sapply(strsplit(s, split = split), cap, USE.NAMES = !is.null(names(s)))
}