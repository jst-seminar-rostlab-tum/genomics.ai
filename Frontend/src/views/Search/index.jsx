import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, Tabs, Tab, CircularProgress } from "@mui/material";
import {
  useHistory,
  useLocation,
  useParams,
  useRouteMatch,
} from "react-router-dom";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchResultList from "components/Search/SearchResultList";
import querySearch from "shared/mock/search";
import TeamCard from "components/Search/SearchResultList/TeamCard";
import InstitutionCard from "components/Search/SearchResultList/InstitutionCard";
import UserCard from "components/Search/SearchResultList/UserCard";
import ProjectCard from "components/Search/SearchResultList/ProjectCard";
import ResultStatus from "components/Search/ResultStatus";
import { Switch } from "react-router-dom";
import { Route } from "react-router-dom";

const Search = ({ sidebarShown }) => {
  /* Booleans */
  const paddingL = useCallback(
    () => (sidebarShown ? "130px" : "380px"),
    [sidebarShown]
  );

  const history = useHistory();
  const location = useLocation();
  const { path } = useRouteMatch();

  // type of the items searched (teams/institutions/users/projects)
  const { type } = useParams();
  const submittedKeyword =
    new URLSearchParams(location.search).get("keyword") || "";

  const [searchedKeyword, setSearchedKeyword] = useState(submittedKeyword);
  const [searchedData, setSearchedData] = useState([]);
  const [isLoading, setIsLoading] = useState(true);

  const searchedKeywordChangeHandler = (event) => {
    setSearchedKeyword(event.target.value);
  };

  const changedTabHandler = (event, newValue) => {
    setIsLoading(true);
    const newPath = setTypeInUrl(path, newValue);
    history.push(`${newPath}`);
  };

  const fetchSearchHandler = useCallback(async (type, keyword) => {
    history.push(`?keyword=${keyword}`);
    const searchResponse = await querySearch(type, keyword.toLowerCase());
    setSearchedData(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    fetchSearchHandler(type, searchedKeyword);
  }, [fetchSearchHandler, type]);

  const submitSearch = () => {
    fetchSearchHandler(type, searchedKeyword);
  };

  const renderSearchResultsList = (listItemWrapper) => {
    return (
      <SearchResultList
        listItemWrapper={listItemWrapper}
        searchedData={searchedData}
      />
    );
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
          <Tabs value={type} onChange={changedTabHandler}>
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
            <React.Fragment>
              <ResultStatus
                count={searchedData.length}
                searchedEntity={type}
                searchedKeyword={submittedKeyword}
              />
              <Switch>
                <Route path={setTypeInUrl(path, "teams")}>
                  {renderSearchResultsList(TeamCard)}
                </Route>
                <Route path={setTypeInUrl(path, "institutions")}>
                  {renderSearchResultsList(InstitutionCard)}
                </Route>
                <Route path={setTypeInUrl(path, "users")}>
                  {renderSearchResultsList(UserCard)}
                </Route>
                <Route path={setTypeInUrl(path, "projects")}>
                  {renderSearchResultsList(ProjectCard)}
                </Route>
              </Switch>
            </React.Fragment>
          )}
        </Box>
      </div>
    </Stack>
  );
};

function createUrl(path, pathParam, value) {
  return path.replace(`:${pathParam}`, value);
}

function setTypeInUrl(path, value) {
  return createUrl(path, "type", value);
}

export default Search;
