import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, CircularProgress } from "@mui/material";
import {
  useHistory,
  useLocation,
  useParams,
  useRouteMatch,
  Route,
} from "react-router-dom";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchTabs from "components/Search/SearchTabs";
import SearchContent from "components/Search/SearchContent";
import Filter from "components/Search/Filter";
import GeneralFilter from "components/Search/Filter/GeneralFilter";
import { setTypeInUrl } from "shared/utils/common/utils";
import querySearch from "shared/mock/search";
import TeamsFilter from "components/Search/Filter/TeamsFilter";

const Search = ({ sidebarShown }) => {
  /* Booleans */
  const paddingL = useCallback(
    () => (sidebarShown ? "130px" : "380px"),
    [sidebarShown]
  );

  // state managed in path and query params
  const history = useHistory();
  const { search } = useLocation();
  const { path } = useRouteMatch();

  const searchParams = new URLSearchParams(search);

  // type of the searched items (teams/institutions/users/projects)
  const { type } = useParams();
  const searchedKeyword = searchParams.get("keyword") || "";

  const [searchedData, setSearchedData] = useState([]);
  const [isLoading, setIsLoading] = useState(true);

  // function to update the state in the URL
  const updateQueryParams = (param, value) => {
    const params = new URLSearchParams(history.location.search);
    if (value) {
      params.set(param, value);
    } else {
      params.delete(param);
    }

    history.push({
      pathname: history.location.pathname,
      search: params.toString(),
    });
  };

  const searchedKeywordChangeHandler = (event) => {
    updateQueryParams("keyword", event.target.value);
  };

  const changedTabHandler = (event, newValue) => {
    setIsLoading(true);
    const newPath = setTypeInUrl(path, newValue);
    history.push({ pathname: `${newPath}`, search: history.location.search });
  };

  const fetchSearchHandler = useCallback(async (type, keyword) => {
    const searchResponse = await querySearch(type, keyword.toLowerCase());
    setSearchedData(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    fetchSearchHandler(type, searchedKeyword);
  }, [fetchSearchHandler, type, searchedKeyword]);

  return (
    <Stack direction="column" sx={{ paddingLeft: paddingL }}>
      <div className={styles.title}>
        <h1>Search</h1>
        <Box sx={{ margin: "auto", maxWidth: 1200 }}>
          <SearchBar
            searchedKeyword={searchedKeyword}
            searchedKeywordChangeHandler={searchedKeywordChangeHandler}
            filterComponent={
              <Filter
                searchParams={searchParams}
                updateQueryParams={updateQueryParams}
                path={path}
              />
            }
          />
          <SearchTabs value={type} onChange={changedTabHandler} />
          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center" }}>
              <CircularProgress />
            </Box>
          )}
          {!isLoading && (
            <SearchContent
              searchedData={searchedData}
              type={type}
              searchedKeyword={searchedKeyword}
            />
          )}
        </Box>
      </div>
    </Stack>
  );
};

export default Search;
