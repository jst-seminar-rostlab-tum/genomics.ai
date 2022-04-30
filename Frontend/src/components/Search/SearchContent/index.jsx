import React from "react";

import { useRouteMatch, Switch, Route } from "react-router-dom";
import SearchResultList from "components/Search/SearchResultList";
import TeamCard from "components/Search/SearchResultList/TeamCard";
import InstitutionCard from "components/Search/SearchResultList/InstitutionCard";
import UserCard from "components/Search/SearchResultList/UserCard";
import ProjectCard from "components/Search/SearchResultList/ProjectCard";
import ResultStatus from "components/Search/ResultStatus";
import { setSeachCategoryInUrl } from "shared/utils/common/utils";


// wrapper component to display the searched items  
const SearchContent = ({ searchResult, searchCategory, searchedKeyword }) => {
  const { path } = useRouteMatch();

  const renderSearchResultsList = (listItemWrapper) => {
    return (
      <SearchResultList
        listItemWrapper={listItemWrapper}
        searchResult={searchResult}
      />
    );
  };

  return (
    <React.Fragment>
      <ResultStatus
        count={searchResult.length}
        searchedEntity={searchCategory}
        searchedKeyword={searchedKeyword}
      />
      <Switch>
        <Route path={setSeachCategoryInUrl(path, "teams")}>
          {renderSearchResultsList(TeamCard)}
        </Route>
        <Route path={setSeachCategoryInUrl(path, "institutions")}>
          {renderSearchResultsList(InstitutionCard)}
        </Route>
        <Route path={setSeachCategoryInUrl(path, "users")}>
          {renderSearchResultsList(UserCard)}
        </Route>
        <Route path={setSeachCategoryInUrl(path, "projects")}>
          {renderSearchResultsList(ProjectCard)}
        </Route>
      </Switch>
    </React.Fragment>
  );
};

export default SearchContent;
