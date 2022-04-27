import React, { useCallback, useState, useEffect } from "react";
import { Box, Stack, CircularProgress } from "@mui/material";
import {
  useHistory,
  useLocation,
  useParams,
  useRouteMatch,
} from "react-router-dom";

import styles from "./search.module.css";
import SearchBar from "components/Search/SearchBar";
import SearchTabs from "components/Search/SearchTabs";
import SearchContent from "components/Search/SearchContent";
import { setTypeInUrl } from "shared/utils/common/utils";
import querySearch from "shared/mock/search";

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
  }, [fetchSearchHandler, type, searchedKeyword]);


  return (
    <Stack direction="column" sx={{ paddingLeft: paddingL }}>
      <div className={styles.title}>
        <h1>Search</h1>
        <Box sx={{ margin: "auto", maxWidth: 1200 }}>
          <SearchBar
            searchedKeyword={searchedKeyword}
            searchedKeywordChangeHandler={searchedKeywordChangeHandler}
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
              submittedKeyword={submittedKeyword}
            />
          )}
        </Box>
      </div>
    </Stack>
  );
};


export default Search;
