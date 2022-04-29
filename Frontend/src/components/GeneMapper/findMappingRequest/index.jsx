import axios from 'axios';
import React, { useState, useEffect } from 'react';
// import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
// import SearchIcon from '@mui/icons-material/Search';
import './findMappingRequest.css';

function FindMappingRequest() {
  const [projects, setProjects] = useState([]);
  const urlDummy = 'https://jsonplaceholder.typicode.com/users';

  //   useEffect(() => {
  // axios.get(`${BACKEND_ADDRESS}/api/projects`).then((res) =>
  //    { setProjects(res.data.results.map((p) => p.name)); });
  //   },
  //   []);

  //   useEffect(() => {
  //     axios.get(urlDummy).then((response) => response.json()).then((json)
  //   => { setProjects(json); });
  //   },
  //   []);
  // fetch datafrom URL
  useEffect(() => {
    const fetchData = async () => {
      const result = await axios(
        urlDummy,
      );

      setProjects(result.data);
    };

    fetchData();
  }, []);

  //   fetch(urlDummy)
  //     .then((response) => response.json())
  //     .then((json) => setProjects(json));

  const PROJECTS = projects;
  console.log(PROJECTS[1]);

  const [foundUsers, setFoundUsers] = useState(PROJECTS);
  const [mappingName, setMappingName] = useState('');
  // console.log(mappingName);

  const filter = (e) => {
    const keyword = e.target.value;

    if (keyword !== '') {
      const results = PROJECTS.filter(
        (owner) => owner.name.toLowerCase().startsWith(keyword.toLowerCase()),
        // Use the toLowerCase() method to make it case-insensitive
      );
      setFoundUsers(results);
    } else {
      setFoundUsers(PROJECTS);
      // If the text field is empty, show all users
    }
    setMappingName(keyword);
  };

  return (
    <div className="formInput">
      <input
        type="search"
        placeholder="Find a mapping"
        name="name of mapping"
        value={mappingName}
        onChange={filter}
      />

      <div className="user-list">
        {foundUsers && foundUsers.length > 0 ? (
          foundUsers.map((user) => (
            <li key={user.id} className="user">
              {/* <span className="user-id">{user.id}</span> */}
              <span className="user-name">{user.name}</span>
              {/* <span className="user-age">
                {user.age}
                {' '}
                year old
              </span> */}
            </li>
          ))
        ) : (
          <h1>No results found!</h1>
        )}
      </div>

    </div>
  );
}

export default FindMappingRequest;
