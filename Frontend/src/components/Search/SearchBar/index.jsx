import React from "react";
import { Button, FormControl, OutlinedInput, InputAdornment, InputLabel } from "@mui/material";


const SearchBar = (props) => {
  return (
    <FormControl sx={{ m: 1, width: "25ch" }} variant="outlined">
      <InputLabel htmlFor="outlined-adornment-search">Type...</InputLabel>
      <OutlinedInput
        id="outlined-adornment-search"
        type={"text"}
        value={props.searchedKeyword}
        onChange={props.searchedKeywordChangeHandler}
        endAdornment={
          <InputAdornment position="end">
            <Button variant="contained" onClick={props.submitSearch}>
              Search
            </Button>
          </InputAdornment>
        }
        label="Search"
      />
    </FormControl>
  );
};

export default SearchBar;
