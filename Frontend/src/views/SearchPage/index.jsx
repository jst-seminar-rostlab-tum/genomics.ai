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
import Search from 'components/Search';
import UserService from 'shared/services/User.service';
import TeamService from 'shared/services/Team.service';
import InstitutionService from 'shared/services/Institution.service';

import { useAuth } from 'shared/context/authContext';

// definitely target to change, when backend will provide full data
async function getTeams(filterParams) {
  const searchResponse = await TeamService.getTeams(filterParams);
  const teamsWithInstitutions = searchResponse.filter((team) => team.institutionId);
  const institutionRequests = teamsWithInstitutions.map(
    (team) => InstitutionService.getInstitutionById(team.institutionId)
    ,
  );
  const institutionsResponse = await Promise.all(institutionRequests);
  institutionsResponse.forEach(
    (institution,
      index) => {
      teamsWithInstitutions[index].institution = institution;
    },
  );
  return searchResponse;
}

// definitely target to change, when backend will provide full data
async function getInstitutions(filterParams) {
  const searchResponse = await InstitutionService.getInstitutions(filterParams);
  const teamsRequests = searchResponse.map(
    (team) => InstitutionService.getTeamsOfInstitutionById(team.id),
  );
  const teamsResponse = await Promise.all(teamsRequests);
  teamsResponse.forEach(
    (team,
      index) => {
      searchResponse[index].teamsCount = team.length;
    },
  );
  return searchResponse;
}

const SearchPage = () => {
  const [user] = useAuth();

  // state managed in path and query params
  const history = useHistory();
  const { search } = useLocation();
  const { path } = useRouteMatch();

  const searchParams = new URLSearchParams(search);

  // category of the searched items (teams/institutions/users)
  const { searchCategory } = useParams();
  const searchedKeyword = searchParams.get('keyword') || '';

  const [searchRequestResult, setSearchRequestResult] = useState([]);
  const [loadedCategory, setLoadedCategory] = useState('');

  // check if searchRequestResult is of the requested category
  const isLoading = loadedCategory !== searchCategory;

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

  const fetchSearchHandler = useCallback(async () => {
    let searchResponse = [];
    const filterParams = Object.fromEntries(new URLSearchParams(search));
    switch (searchCategory) {
      case 'users':
        searchResponse = await UserService.getUsers(filterParams);
        break;
      case 'teams':
        searchResponse = await getTeams(filterParams);
        break;
      case 'institutions':
        searchResponse = await getInstitutions(filterParams);
        break;
      default:
    }
    setSearchRequestResult(searchResponse);
    setLoadedCategory(searchCategory);
  }, [searchCategory, search]);

  useEffect(() => {
    fetchSearchHandler();
  }, [fetchSearchHandler]);

  return (
    <Stack direction="column" sx={{ paddingLeft: '130px' }}>
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
            value={searchedKeyword}
          />
          <SearchTabs
            value={searchCategory}
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
              user={user}
              fetchSearchHandler={fetchSearchHandler}
            />
          )}
        </Box>
      </div>
    </Stack>
  );
};

export default SearchPage;
