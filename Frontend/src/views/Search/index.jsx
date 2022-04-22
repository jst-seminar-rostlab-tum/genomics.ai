import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, Tabs, Tab, CircularProgress } from "@mui/material";

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

  const [selectedTab, setSelectedTab] = useState("Teams");
  const [searchedKeyword, setSearchedKeyword] = useState("");
  const [searchedData, setSearchedData] = useState([]);
  const [isLoading, setIsLoading] = useState(true);

  const searchedKeywordChangeHandler = (event) => {
    setSearchedKeyword(event.target.value);
  };

  const changedTabHandler = (event, newValue) => {
    setSelectedTab(newValue);
    setIsLoading(true);
  };

  const submitSearch = () => {
    console.log(searchedKeyword);
  };

  const fetchSearchHandler = useCallback(async (type) => {
    const searchResponse = await querySearch(type);
    setSearchedData(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    fetchSearchHandler(selectedTab);
  }, [fetchSearchHandler, selectedTab]);

  console.log(isLoading);
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
          <Tabs value={selectedTab} onChange={changedTabHandler}>
            <Tab label="Teams" value="Teams" />
            <Tab label="Institutions" value="Institutions" />
          </Tabs>
          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center" }}>
              <CircularProgress />
            </Box>
          )}
          {!isLoading && (
            <SearchResultList searchedData={searchedData} type={selectedTab} />
          )}
        </Box>
      </div>
    </Stack>
  );
};

export default Search;
