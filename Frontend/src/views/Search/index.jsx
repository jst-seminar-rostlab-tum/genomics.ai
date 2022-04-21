import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack } from "@mui/material";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchResultList from "components/Search/SearchResultList";
import querySearch from "shared/mock/search";

const Search = ({ sidebarShown }) => {
  /* Booleans */
  const paddingL = useCallback(
    () => (sidebarShown ? "130px" : "380px"),
    [sidebarShown]
  );

  const [searchedKeyword, setSearchedKeyword] = useState("");
  const [searchedData, setSearchedData] = useState([]);

  const searchedKeywordChangeHandler = (event) => {
    setSearchedKeyword(event.target.value);
  };
  const submitSearch = () => {
    console.log(searchedKeyword);
  };

  const fetchSearchHandler = useCallback(async () => {
    const searchResponse = await querySearch();
    setSearchedData(searchResponse);
  }, []);

  useEffect(() => {
    fetchSearchHandler();
  }, [fetchSearchHandler]);

  return (
    <Stack direction="column" sx={{ paddingLeft: paddingL }}>
      {" "}
      <div className={styles.title}>
        <h1>Search</h1>
        <Box sx={{ margin: "auto", maxWidth: 1200 }}>
          <SearchBar
            searchedKeyword={searchedKeyword}
            searchedKeywordChangeHandler={searchedKeywordChangeHandler}
            submitSearch={submitSearch}
          />
          <SearchResultList searchedData={searchedData} />
        </Box>
      </div>
    </Stack>
  );
};

export default Search;
