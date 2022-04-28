import React from "react";
import { FormControl, OutlinedInput, InputLabel, Stack } from "@mui/material";

const SearchBar = (props) => {
  return (
    <Stack direction="row">
      <FormControl sx={{ m: 1, width: "25ch" }} variant="outlined">
        <InputLabel htmlFor="outlined-adornment-search">Type...</InputLabel>
        <OutlinedInput
          id="outlined-adornment-search"
          type={"text"}
          value={props.searchedKeyword}
          onChange={props.searchedKeywordChangeHandler}
          label="Search"
        />
      </FormControl>
      {props.filterComponent}
    </Stack>
  );
};

export default SearchBar;
