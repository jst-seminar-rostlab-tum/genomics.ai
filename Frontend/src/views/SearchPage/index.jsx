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
import ProjectService from 'shared/services/Project.service';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';

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
      teamsWithInstitutions[index].institutionTitle = institution.name;
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

  const fetchSearchHandler = useCallback(async (_searchCategory, _searchParams) => {
    let searchResponse = [];
    const filterParams = Object.fromEntries(new URLSearchParams(_searchParams));
    switch (_searchCategory) {
      case 'users':
        searchResponse = await UserService.getUsers(filterParams);
        break;
      case 'teams':
        searchResponse = await getTeams(filterParams);
        break;
      case 'institutions':
        searchResponse = await getInstitutions(filterParams);
        break;
      case 'projects':
        searchResponse = await ProjectService.getProjects(filterParams);
        break;
      case 'atlases':
        searchResponse = await AtlasService.getAtlases();
        break;
      case 'models':
        searchResponse = await ModelService.getModels();
        break;
      default:
    }
    setSearchRequestResult(searchResponse);
    setIsLoading(false);
  }, []);

  useEffect(() => {
    setIsLoading(true);
    fetchSearchHandler(searchCategory, search);
  }, [fetchSearchHandler, searchCategory, search]);

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
