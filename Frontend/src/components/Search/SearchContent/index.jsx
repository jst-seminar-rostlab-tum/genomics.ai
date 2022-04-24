import React from "react";

import { useRouteMatch, Switch, Route } from "react-router-dom";
import SearchResultList from "components/Search/SearchResultList";
import TeamCard from "components/Search/SearchResultList/TeamCard";
import InstitutionCard from "components/Search/SearchResultList/InstitutionCard";
import UserCard from "components/Search/SearchResultList/UserCard";
import ProjectCard from "components/Search/SearchResultList/ProjectCard";
import ResultStatus from "components/Search/ResultStatus";
import { setTypeInUrl } from "shared/utils/common/utils";



const SearchContent = ({ searchedData, type, submittedKeyword }) => {
  const { path } = useRouteMatch();

  const renderSearchResultsList = (listItemWrapper) => {
    return (
      <SearchResultList
        listItemWrapper={listItemWrapper}
        searchedData={searchedData}
      />
    );
  };

  return (
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
  );
};

export default SearchContent;
