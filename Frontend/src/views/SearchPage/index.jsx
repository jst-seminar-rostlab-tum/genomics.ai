import React, { useCallback, useState, useEffect } from 'react';
import { Box, Stack, CircularProgress } from '@mui/material';
import {
  useHistory,
  useLocation,
  useParams,
  useRouteMatch,
} from 'react-router-dom';

import styles from './search.module.css';
import SearchTabs from 'components/SearchPageComponents/SearchTabs';
import SearchContent from 'components/SearchPageComponents/SearchContent';
import Filter from 'components/SearchPageComponents/Filter';
// import { setSeachCategoryInUrl } from "shared/utils/common/utils";
import querySearch from 'shared/mock/search';
import Search from 'components/Search';

const SearchPage = ({ sidebarShown }) => {
  /* Booleans */
  const paddingL = () => (sidebarShown ? '130px' : '380px');

  // state managed in path and query params
  const history = useHistory();
  const { search } = useLocation();
  const { path } = useRouteMatch();

  const searchParams = new URLSearchParams(search);

  // category of the searched items (teams/institutions/users/projects)
  const { searchCategory } = useParams();
  const searchedKeyword = searchParams.get('keyword') || '';

  const [searchRequestResult, setSearchRequestResult] = useState([]);
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

  const searchedKeywordChangeHandler = (value) => {
    updateQueryParams('keyword', value);
  };

  const changedTabHandler = () => {
    setIsLoading(true);
  };

  const fetchSearchHandler = useCallback(async (_searchCategory, keyword) => {
    const searchResponse = await querySearch(
      _searchCategory,
      keyword.toLowerCase(),
    );
    setSearchRequestResult(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    setIsLoading(true);
    fetchSearchHandler(searchCategory, searchedKeyword);
  }, [fetchSearchHandler, searchCategory, searchedKeyword]);

  return (
    <Stack direction="column" sx={{ paddingLeft: paddingL }}>
      <div className={styles.title}>
        <h1>Search</h1>
        <Box sx={{ margin: 'auto', maxWidth: 1200 }}>
          <Search
            filterComponent={(
              <Filter
                searchParams={searchParams}
                updateQueryParams={updateQueryParams}
                path={path}
              />
            )}
            handleSearch={searchedKeywordChangeHandler}
            value={searchedKeyword} // currently two-way-binding missing
          />
          <SearchTabs
            value={searchCategory}
            onChange={changedTabHandler}
            searchParams={searchParams}
            path={path}
          />
          {isLoading && (
            <Box sx={{ display: 'flex', justifyContent: 'center' }}>
              <CircularProgress />
            </Box>
          )}
          {!isLoading && (
            <SearchContent
              searchResult={searchRequestResult}
              searchCategory={searchCategory}
              searchedKeyword={searchedKeyword}
            />
          )}
        </Box>
      </div>
    </Stack>
  );
};

export default SearchPage;
