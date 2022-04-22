import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, Tabs, Tab, CircularProgress } from "@mui/material";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchResultList from "components/Search/SearchResultList";
import querySearch from "shared/mock/search";
import TeamCard from "components/Search/SearchResultList/TeamCard";
import InstitutionCard from "components/Search/SearchResultList/InstitutionCard";

const Search = ({ sidebarShown }) => {
  /* Booleans */
  const paddingL = useCallback(
    () => (sidebarShown ? "130px" : "380px"),
    [sidebarShown]
  );

  const [selectedTab, setSelectedTab] = useState("teams");
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

  const listItemWrapper = {
    teams: TeamCard,
    institutions: InstitutionCard,
  };

  return (
    <Stack direction="column" sx={{ paddingLeft: paddingL }}>
      <div className={styles.title}>
        <h1>Search</h1>
        <Box sx={{ margin: "auto", maxWidth: 1200 }}>
          <SearchBar
            searchedKeyword={searchedKeyword}
            searchedKeywordChangeHandler={searchedKeywordChangeHandler}
            submitSearch={submitSearch}
          />
          <Tabs value={selectedTab} onChange={changedTabHandler}>
            <Tab label="Teams" value="teams" />
            <Tab label="Institutions" value="institutions" />
          </Tabs>
          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center" }}>
              <CircularProgress />
            </Box>
          )}
          {!isLoading && (
            <SearchResultList
              listItemWrapper={listItemWrapper[selectedTab]}
              searchedData={searchedData}
            />
          )}
        </Box>
      </div>
    </Stack>
  );
};

export default Search;
