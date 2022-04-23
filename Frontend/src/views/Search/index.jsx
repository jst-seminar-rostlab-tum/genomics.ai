import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, Tabs, Tab, CircularProgress } from "@mui/material";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchResultList from "components/Search/SearchResultList";
import querySearch from "shared/mock/search";
import TeamCard from "components/Search/SearchResultList/TeamCard";
import InstitutionCard from "components/Search/SearchResultList/InstitutionCard";
import UserCard from "components/Search/SearchResultList/UserCard";
import ProjectCard from "components/Search/SearchResultList/ProjectCard";

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

  const fetchSearchHandler = useCallback(async (type, keyword) => {
    const searchResponse = await querySearch(type, keyword);
    setSearchedData(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    fetchSearchHandler(selectedTab, "");
  }, [fetchSearchHandler, selectedTab]);

  const submitSearch = () => {
    fetchSearchHandler(selectedTab, searchedKeyword);
  };

  const listItemWrapper = {
    teams: TeamCard,
    institutions: InstitutionCard,
    users: UserCard,
    projects: ProjectCard,
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
            <Tab label="Users" value="users" />
            <Tab label="Projects" value="projects" />
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
